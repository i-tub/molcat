"""
Display a 2D sketch of a structure, from a SMILES or a file, to a terminal that
supports graphics, such as kitty, Ghostty, and iTerm2.
"""

__version__ = '0.2.0'

import array
import argparse
import base64
import fcntl
import gzip
import itertools
import logging
import os
import re
import sys
import termios
import tty
from pathlib import Path
from typing import Any, Iterator, Iterable

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdkit_logger = logging.getLogger('rdkit')

MAX_ATOMS = 500

# Y/X shape of the generated image, unless -size_y is specified.
ASPECT_RATIO = 3 / 5

# By default, use image width of DEFAULT_REL_X * window_width
DEFAULT_REL_X = 0.8

# Size to use by default if we can't determine window width.
DEFAULT_SIZE = (500, 300)

MolSupplier = Chem.SmilesMolSupplier | Chem.SDMolSupplier | Chem.MaeMolSupplier


def get_window_size():
    """
    Return the window size (width, height) in pixels.
    """
    buf = array.array('H', [0, 0, 0, 0])
    fcntl.ioctl(sys.stdout, termios.TIOCGWINSZ, buf)
    _, _, width, height = buf
    return width, height


def _chunk(data, chunk_size=2048):
    size = len(data)
    for start in range(0, size, chunk_size):
        chunk = data[start:start + chunk_size]
        last = start + chunk_size >= size
        yield chunk, last


def _show_image_chunked(png_data: bytes) -> None:
    """
    Print a PNG file to the terminal using the Kitty protocol.
    """
    b64_data = base64.b64encode(png_data)
    # a=T (Transfer & Display), f=100 (PNG), q=2 (No Acks/Gibberish)
    cmd = b'a=T,f=100,q=2,m='
    for chunk, last in _chunk(b64_data):
        flag = b'0' if last else b'1'
        sys.stdout.buffer.write(b'\033_G' + cmd + flag + b';' + chunk +
                                b'\033\\')
    print(flush=True)


def show_image(png_data: bytes) -> None:
    """
    Print a PNG file to the terminal using the Kitty protocol.
    """
    if 'kitty' in os.environ.get('TERM', ''):
        return _show_image_chunked(png_data)
    b64_data = base64.b64encode(png_data)
    # a=T (Transfer & Display), f=100 (PNG), q=2 (No Acks/Gibberish)
    cmd = b'a=T,f=100,q=2'
    sys.stdout.buffer.write(b'\033_G' + cmd + b';' + b64_data + b'\033\\')
    print(flush=True)


def get_png(mol: Chem.Mol, size: tuple[int, int] = DEFAULT_SIZE) -> bytes:
    """
    Draw a molecule to the terminal.
    """
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts: Any = d.drawOptions()
    opts.addStereoAnnotation = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def get_thumbnails_png(mols,
                       cols: int,
                       thumbnail_size: tuple[int, int],
                       legend_prop: str = 'index'):

    mols = list(mols)
    if legend_prop != 'none':
        legends = [
            mol.GetProp(legend_prop) if mol.HasProp(legend_prop) else ''
            for mol in mols
        ]
    else:
        legends = None
    return Draw.MolsToGridImage(mols,
                                returnPNG=True,
                                molsPerRow=cols,
                                subImgSize=thumbnail_size,
                                legends=legends)


def show_mol(mol: Chem.Mol, size: tuple[int, int] = DEFAULT_SIZE) -> None:
    """
    Draw a molecule to the terminal.
    """
    show_image(get_png(mol, size))


def copy_mol(mol: Chem.Mol, size: tuple[int, int] = DEFAULT_SIZE) -> None:
    """
    Copy molecule image to the clipboard using kitty protocol.
    """
    copy_5522(get_png(mol, size), 'image/png')


def copy_5522(data: bytes, mime_type: str) -> None:
    """
    Copy a binary object to the clipboard using kitty protocol (OSC 5522).
    """
    OSC5522 = b'\033]5522;'
    ST = b'\033\\'
    b64_data = base64.b64encode(data)
    b64_type = base64.b64encode(mime_type.encode('ascii'))

    sys.stdout.buffer.write(OSC5522 + b'type=write' + ST)
    for chunk, _ in _chunk(b64_data):
        sys.stdout.buffer.write(OSC5522 + b'type=wdata:mime=' + b64_type +
                                b';' + chunk + ST)
    sys.stdout.buffer.write(OSC5522 + b'type=wdata' + ST)
    sys.stdout.flush()

    try:
        with open('/dev/tty', 'rb+', buffering=0) as tty_file:
            fd = tty_file.fileno()
            old_settings = termios.tcgetattr(fd)
            try:
                # Set to raw mode so we don't wait for a newline from the terminal
                tty.setraw(fd)
                response = b""
                while True:
                    char = tty_file.read(1)
                    response += char
                    # Kitty's response ends with the String Terminator
                    if response.endswith(ST) or response.endswith(b'\a'):
                        break
                # Optional: Check 'status=DONE' in response for error handling
            finally:
                termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    except (FileNotFoundError, OSError):
        # Fallback for environments without a TTY (like CI/CD or cron)
        pass


def to_2d(mol: Chem.Mol,
          keeph: bool = False,
          idx: int | None = None,
          cleanIt: bool = True) -> Chem.Mol:
    """
    Return a modified copy of mol with 2D coordinates added, and depending on
    flags, hydrogens removed and indices added as atomNote properties. The notes
    have the original atom indices, before deleting the hydrogens.

    :param keeph: if False, remove the hydrogen atoms
    :param idx: if not None, set the atomNote property of each atom to its atom
                index plus the value of `idx`. This is normally 0 or 1,
                depending on which numbering convention is desired.
    """
    mol = Chem.Mol(mol)
    if idx is not None:
        for atom in mol.GetAtoms():
            atom.SetProp('atomNote', str(atom.GetIdx() + idx))
    if not keeph:
        mol = Chem.RemoveHs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=cleanIt, force=True)

    rdDepictor.Compute2DCoords(mol)
    Chem.Kekulize(mol, True)
    return mol


def get_reader(filename, removeHs=False) -> MolSupplier:
    """
    Return a Mol supplier for the given filename.
    """
    if filename.endswith('.smi'):
        return Chem.SmilesMolSupplier(filename, titleLine=False)
    elif filename.endswith('.csv'):
        return Chem.SmilesMolSupplier(filename, delimiter=',')
    elif filename.endswith('.sdf') or filename.endswith('.mol'):
        return Chem.SDMolSupplier(filename, removeHs=removeHs)
    elif filename.endswith('.mae'):
        return Chem.MaeMolSupplier(filename, removeHs=removeHs)
    elif filename.endswith('.maegz') or filename.endswith('.mae.gz'):
        return Chem.MaeMolSupplier(gzip.open(filename), removeHs=removeHs)
    else:
        raise ValueError(f'Unknown file format for {filename}')


def parse_range(n: str) -> tuple[int, int]:
    """
    Parse a range string formated as "<start>-<stop>", but note that the range
    string follows a start from 1 and an end-inclusive convention, so "1-3"
    would become `(0, 3)`.
    """
    try:
        start = int(n)
        return (start - 1, start)
    except ValueError:
        start, stop = (int(s) for s in n.split('-'))
        return start - 1, stop


def LogLevel(level_name: str) -> int:
    try:
        return getattr(logging, level_name)
    except AttributeError:
        raise TypeError(f'Invalid log level: {level_name}')


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file_or_smiles',
                        nargs='?',
                        help="structure input file or SMILES strings. "
                        "If not provided, SMILES will be read from stdin.")
    parser.add_argument(
        '-n',
        default='1',
        help="index of structure to display. May be a range ('-n 1-4'). "
        "Default: 1 (only show the first structure in the file)")
    parser.add_argument('--all',
                        '-a',
                        action='store_true',
                        help='show all structures in the file')

    # args.idx is implemented in a slightly hacky way: two different flags
    # have it as their destination: one sets it to 1, the other to 0, and the
    # default is None.
    idx_group = parser.add_mutually_exclusive_group()
    idx_group.add_argument('--idx',
                           '-i',
                           action='store_true',
                           default=None,
                           help='show atom indexes (1-based)')
    idx_group.add_argument('--zidx',
                           '-z',
                           dest='idx',
                           action='store_const',
                           const=0,
                           default=None,
                           help='show atom indexes (0-based)')

    parser.add_argument('--keeph',
                        '-H',
                        action='store_true',
                        help='keep all hydrogen atoms')
    parser.add_argument(
        '--size-x',
        '-x',
        type=int,
        help='X dimension in pixels; default: determine automatically')
    parser.add_argument(
        '--size-y',
        '-y',
        type=int,
        default=None,
        help='X dimension in pixels; default is a function of -x')
    parser.add_argument('--cols', type=int, help='number of columns to use')
    parser.add_argument(
        '--legend',
        default='index',
        help='legend to use when using --cols. Default="index". Mol name is '
        'shown with "name" or "title". Anything else is a property.')
    parser.add_argument(
        '--min-atoms',
        type=int,
        default=1,
        help='only display molecules with at least this many atoms')
    parser.add_argument('--log-level',
                        type=LogLevel,
                        default=logging.FATAL,
                        help='RDKit log level; default="FATAL"')
    parser.add_argument('--copy',
                        '-c',
                        action='store_true',
                        help='copy to clipboard')
    parser.add_argument('--write',
                        '-w',
                        metavar='<file.png>',
                        help='Write to PNG file')
    args = parser.parse_args(argv)
    if args.legend in ('name', 'title'):
        args.legend = '_Name'

    return args


def mols_from_str(line: str,
                  strict: bool = False,
                  min_atoms: int = 1) -> Iterator[Chem.Mol]:
    """
    Split a line on comma or whitespace into SMILES strings, ignoring any that
    are not valid.
    """
    for tok in re.split(r'''[,\s"']+''', line):
        # We don't sanitize during the initial conversion so we can keep any
        # graph hydrogens until to_2d decides whether to remove them or not.
        mol = Chem.MolFromSmiles(tok, sanitize=False)
        if mol is not None:
            if mol.GetNumAtoms() >= min_atoms:
                Chem.SanitizeMol(mol)
                yield mol
        elif strict:
            raise ValueError(f'Invalid SMILES: {tok}')


def mols_from_file(file: Iterable[str],
                   strict: bool = False,
                   min_atoms: int = 1) -> Iterator[Chem.Mol]:
    for line in file:
        yield from mols_from_str(line.strip(), strict, min_atoms)


def get_mols(args):
    if args.file_or_smiles and os.path.isfile(args.file_or_smiles):
        reader = get_reader(args.file_or_smiles)
        if args.all:
            start = 0
        else:
            start, stop = parse_range(args.n)
            reader = itertools.islice(reader, start, stop)
        for i, mol in enumerate(reader, start + 1):
            if not mol or mol.GetNumAtoms() > MAX_ATOMS:
                continue
            mol.SetIntProp('index', i)
            yield mol
    else:
        # SMILES text mode. This one does not support specifying a range,
        # because it allows for some sloppiness, just using any "words" that
        # from stdin or the command line that happen to be valid SMILES.
        strict = bool(args.file_or_smiles)
        lines = [args.file_or_smiles] if args.file_or_smiles else sys.stdin
        yield from mols_from_file(lines, strict, args.min_atoms)


def determine_size(x: int = 0, y: int = 0):
    """
    Determine size for images. If only one of x or y is provided, apply the
    default aspect ratio. If neither is provided, determine the size
    automatically based on the window width.
    """
    if x and y:
        return x, y
    elif x and not y:
        return x, int(x * ASPECT_RATIO)
    elif y and not x:
        return int(y / ASPECT_RATIO), y
    else:
        x, _ = get_window_size()
        if x:
            x = int(x * DEFAULT_REL_X)
            return int(x), int(x * ASPECT_RATIO)
        else:
            return DEFAULT_SIZE


def main():
    args = parse_args()

    # Capture RDKit warnings
    rdkit_logger.setLevel(args.log_level)
    Chem.rdBase.LogToPythonLogger()

    try:
        mols = get_mols(args)
        mols2d = (to_2d(mol, args.keeph, args.idx) for mol in mols)
        if args.cols:
            x = args.size_x or get_window_size()[0] or DEFAULT_SIZE[0]
            thumbnail_size = (int(x / args.cols),
                              int(x / args.cols * ASPECT_RATIO))
            png_data = get_thumbnails_png(mols2d, args.cols, thumbnail_size,
                                          args.legend)
            show_image(png_data)
        else:
            size = determine_size(args.size_x, args.size_y)
            for mol in mols2d:
                png_data = get_png(mol, size)
                show_image(png_data)
        if args.copy:
            copy_5522(png_data, 'image/png')
        if args.write:
            Path(args.write).write_bytes(png_data)
    except ValueError as e:
        sys.exit(e)
