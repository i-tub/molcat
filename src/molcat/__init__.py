"""
Display a 2D sketch of a structure from SMILES or a file to a terminal that
support graphics, such as kitty, iterm, or ghostty.
"""

__version__ = '0.1.0'

import argparse
import base64
import gzip
import os
import sys

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

MAX_ATOMS = 500


def show_image(png_data):
    b64_data = base64.b64encode(png_data)
    # a=T (Transfer & Display), f=100 (PNG), q=2 (No Acks/Gibberish)
    cmd = b'a=T,f=100,q=2'
    sys.stdout.buffer.write(b'\033_G' + cmd + b';' + b64_data + b'\033\\')
    print(flush=True)


def show_structure(mol, size=300, indexes=False):
    d = rdMolDraw2D.MolDraw2DCairo(size, size)
    opts = d.drawOptions()
    d.DrawMolecule(mol)
    d.FinishDrawing()
    show_image(d.GetDrawingText())


def to_2d(mol: Chem.Mol,
          keeph: bool = False,
          idx: int | None = None) -> Chem.Mol:
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
        for bond in mol.GetBonds():
            bond.SetBondDir(Chem.BondDir.NONE)
        Chem.AssignStereochemistryFrom3D(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    rdDepictor.Compute2DCoords(mol)
    return mol


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('file_or_smiles',
                        help="structure input file or SMILES string")
    parser.add_argument(
        '-n',
        default=1,
        help="index of structure to display. May be a range ('-n 1-4')")
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
    parser.add_argument('-size', type=int, default=500)
    return parser.parse_args()


def parse_range(n):
    try:
        start = int(n)
        return start, start
    except ValueError:
        return (int(s) for s in n.split('-'))


def get_reader(filename, index=1):
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
        sys.exit(f'Unknown file format for {filename}')


def main():
    args = parse_args()
    if os.path.isfile(args.file_or_smiles):
        start, end = parse_range(args.n)
        reader = get_reader(args.file_or_smiles)
        for i, mol in enumerate(reader, 1):
            if not mol:
                continue
            if i < start:
                continue
            if i > end and not args.all:
                break
            if mol.GetNumAtoms() > MAX_ATOMS:
                continue
            mol = to_2d(mol, args.keeph, args.idx)
            show_structure(mol, args.size, args.idx)
    else:
        mol = Chem.MolFromSmiles(args.file_or_smiles, sanitize=False)
        Chem.SanitizeMol(mol)
        mol = to_2d(mol, args.keeph, args.idx)
        show_structure(mol, args.size, args.idx)
