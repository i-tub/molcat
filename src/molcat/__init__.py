"""
Display a 2D sketch of a structure, from a SMILES or a file, to a terminal that
support graphics, such as kitty, Ghostty, and iTerm2.
"""

__version__ = '0.1.1'

import argparse
import base64
import gzip
import itertools
import logging
import os
import re
import sys
from typing import Any, Iterator, Iterable

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

rdkit_logger = logging.getLogger('rdkit')

MAX_ATOMS = 500

# Y/X shape of the generated image, unless -size_y is specified.
ASPECT_RATIO = 3 / 5

MolSupplier = Chem.SmilesMolSupplier | Chem.SDMolSupplier | Chem.MaeMolSupplier


def show_image(png_data: bytes) -> None:
    """
    Print a PNG file to the terminal using the Kitty protocol.
    """
    b64_data = base64.b64encode(png_data)
    # a=T (Transfer & Display), f=100 (PNG), q=2 (No Acks/Gibberish)
    cmd = b'a=T,f=100,q=2'
    sys.stdout.buffer.write(b'\033_G' + cmd + b';' + b64_data + b'\033\\')
    print(flush=True)


def show_mol(mol: Chem.Mol, size: tuple[int, int] = (500, 300)) -> None:
    """
    Draw a molecule to the terminal.
    """
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts: Any = d.drawOptions()
    opts.addStereoAnnotation = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    show_image(d.GetDrawingText())


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


def get_reader(filename) -> MolSupplier:
    """
    Return a Mol supplier for the given filename.
    """
    if filename.endswith('.smi'):
        return Chem.SmilesMolSupplier(filename, titleLine=False)
    elif filename.endswith('.csv'):
        return Chem.SmilesMolSupplier(filename, delimiter=',')
    elif filename.endswith('.sdf') or filename.endswith('.mol'):
        return Chem.SDMolSupplier(filename, removeHs=False)
    elif filename.endswith('.mae'):
        return Chem.MaeMolSupplier(filename, removeHs=False)
    elif filename.endswith('.maegz') or filename.endswith('.mae.gz'):
        return Chem.MaeMolSupplier(gzip.open(filename), removeHs=False)
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
    parser.add_argument('--size-x',
                        '-x',
                        type=int,
                        default=500,
                        help='X dimension in pixels; default=%(default)s')
    parser.add_argument(
        '--size-y',
        '-y',
        type=int,
        default=None,
        help='X dimension in pixels; default is a function of -x')
    parser.add_argument(
        '--min-atoms',
        type=int,
        default=1,
        help='only display molecules with at least this many atoms')
    parser.add_argument('--log-level',
                        type=LogLevel,
                        default=logging.FATAL,
                        help='RDKit log level; default="FATAL"')
    args = parser.parse_args(argv)
    args.size_y = args.size_y or int(args.size_x * ASPECT_RATIO)

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
        if not args.all:
            start, stop = parse_range(args.n)
            reader = itertools.islice(reader, start, stop)
        for mol in reader:
            if not mol or mol.GetNumAtoms() > MAX_ATOMS:
                continue
            yield mol
    else:
        # SMILES text mode. This one does not support specifying a range, because
        # it allows for some sloppiness, just using any "words" that from stdin
        # or the command line that happen to be valid SMILES.
        strict = bool(args.file_or_smiles)
        lines = [args.file_or_smiles] if args.file_or_smiles else sys.stdin
        yield from mols_from_file(lines, strict, args.min_atoms)


def main():
    args = parse_args()

    # Capture RDKit warnings
    rdkit_logger.setLevel(args.log_level)
    Chem.rdBase.LogToPythonLogger()

    try:
        for mol in get_mols(args):
            mol2d = to_2d(mol, args.keeph, args.idx)
            show_mol(mol2d, (args.size_x, args.size_y))
    except ValueError as e:
        sys.exit(e)
