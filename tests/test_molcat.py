import molcat

import pytest
from rdkit import Chem


@pytest.fixture()
def args():
    return molcat.parse_args([])


def test_mols_from_file():
    lines = ['hello OCCO world c1ccccc1', 'N the end']

    mols = list(molcat.mols_from_file(lines))

    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['OCCO', 'c1ccccc1', 'N']


def test_mols_from_file_strict():
    with pytest.raises(ValueError, match='Invalid SMILES'):
        list(molcat.mols_from_file(['O', 'x'], strict=True))


def test_get_mols_from_smiles(args):
    args.file_or_smiles = 'O'
    mols = list(molcat.get_mols(args))
    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['O']


def test_get_mols_from_file_default(args):
    # By default, we should only get the first structure in the file.
    args.file_or_smiles = 'tests/test.sdf'
    mols = list(molcat.get_mols(args))
    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['C']


def test_get_mols_from_file_specific_structure(args):
    # We'll ask for structure 2. The file has the first 5 linear alkanes.
    args.file_or_smiles = 'tests/test.sdf'
    args.n = '2'
    mols = list(molcat.get_mols(args))
    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['CC']


def test_get_mols_from_file_range(args):
    # We'll ask for structures 2-3. The file has the first 5 linear alkanes.
    args.file_or_smiles = 'tests/test.sdf'
    args.n = '2-3'
    mols = list(molcat.get_mols(args))
    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['CC', 'CCC']


def test_get_mols_from_file_all(args):
    # We'll ask for structures 2-3. The file has the first 5 linear alkanes.
    args.file_or_smiles = 'tests/test.sdf'
    args.all = True
    mols = list(molcat.get_mols(args))
    got_smiles = [Chem.MolToSmiles(mol) for mol in mols]
    assert got_smiles == ['C', 'CC', 'CCC', 'CCCC', 'CCCCC']
