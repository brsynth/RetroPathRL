#!/usr/bin/env python
import pytest

from utilities.chemtools.Filters import Filters
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import MolFromInchi, MolToInchi


def test_init():
    assert Filters()

def test_copy_properties():
    # TODO: add some tests here
    pass

def test_keep_biggest():
    mol = Filters.keep_biggest(MolFromSmiles('CCCC.CC'))
    assert MolToSmiles(mol) == 'CCCC'
    mol = Filters.keep_biggest(MolFromSmiles('CCCCC.CC.[H].CCC'))
    assert MolToSmiles(mol) == 'CCCCC'
    mol = Filters.keep_biggest(MolFromInchi('InChI=1S/C5H12N2O2.C4H7NO4/c6-3-1-2-4(7)5(8)9;5-2(4(8)9)1-3(6)7/h4H,1-3,6-7H2,(H,8,9);2H,1,5H2,(H,6,7)(H,8,9)/t4-;2-/m00/s1'))
    assert MolToInchi(mol) == 'InChI=1S/C4H7NO4/c5-2(4(8)9)1-3(6)7/h2H,1,5H2,(H,6,7)(H,8,9)/t2-/m0/s1'
    mol = Filters.keep_biggest(MolFromInchi('InChI=1S/Mo.4O/q;;;2*-1'))
    assert MolToInchi(mol) == 'InChI=1S/Mo'

def test_commute_inchi():
    inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1'
    mol = Filters.commute_inchi(MolFromInchi(inchi))
    assert MolToInchi(mol) == inchi
    
def test_remove_isotope():
    mol = Filters.remove_isotope(MolFromSmiles('c1cc[14cH]cc1'))
    assert MolToSmiles(mol) == ('c1ccccc1')

def test_neutralise_charge():
    mol = Filters.neutralise_charge(MolFromSmiles('CC(C(=O)[O-])O'))
    assert MolToSmiles(mol) == ('CC(O)C(=O)O')

def test_add_hydrogen():
    mol = Filters.add_hydrogen(MolFromSmiles('CC(O)C(=O)O'))
    assert MolToSmiles(mol) == '[H]OC(=O)C([H])(O[H])C([H])([H])[H]'
    mol = Filters.add_hydrogen(MolFromSmiles('CC(C(=O)[O-])O'))
    assert MolToSmiles(mol) == '[H]OC([H])(C(=O)[O-])C([H])([H])[H]'
    
def test_kekulize():
    mol = Filters.kekulize(MolFromSmiles('c1ccccc1'))
    assert MolToSmiles(mol) == 'C1=CC=CC=C1'

def test_remove_stereo():
    mol = Filters.remove_stereo(MolFromSmiles('C[C@@H](C(=O)[O-])O'))
    assert MolToSmiles(mol) == 'CC(O)C(=O)[O-]'
    mol = Filters.remove_stereo(MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'))
    assert MolToSmiles(mol) == 'OC1=NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(O)=Nc2ccccc21'
    mol = Filters.commute_inchi(mol)  # Expected to change tautomerism
    assert MolToSmiles(mol) == 'O=C1NC(C2=CNC3=C2C=C(O)C=C3)=CC1=C1C(=O)NC2=CC=CC=C21'
