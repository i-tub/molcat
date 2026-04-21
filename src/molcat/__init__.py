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


def to_2d(mol, keeph=False):
    if keeph:
        mol = Chem.AddHs(mol)
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
    parser.add_argument('-all',
                        action='store_true',
                        help='show all structures in the file')
    parser.add_argument('-idx', action='store_true', help='show atom indexes')
    parser.add_argument('-zidx', action='store_true', help='show atom indexes')
    parser.add_argument('-keeph',
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
        return Chem.SDMolSupplier(filename)
    elif filename.endswith('.mae'):
        return Chem.MaeMolSupplier(filename)
    elif filename.endswith('.maegz') or filename.endswith('.mae.gz'):
        return Chem.MaeMolSupplier(gzip.open(filename))
    else:
        sys.exit(f'Unknown file format for {filename}')


def add_index_props(mol, from_one: bool):
    for atom in mol.GetAtoms():
        atom.SetProp('atomNote', str(atom.GetIdx() + from_one))


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
            if args.idx or args.zidx:
                add_index_props(mol, args.idx)
            mol = to_2d(mol, args.keeph)
            show_structure(mol, args.size, args.idx)
    else:
        mol = Chem.MolFromSmiles(args.file_or_smiles)
        rdDepictor.Compute2DCoords(mol)
        if args.idx or args.zidx:
            add_index_props(mol, args.idx)
        show_structure(mol, args.size, args.idx)
