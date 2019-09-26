#!/usr/bin/env python
import pytest

import inspect
from utilities.chemtools.Standardizer import Standardizer
from utilities.chemtools.Sequences import sequence_tunable
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem import MolFromInchi, MolToInchi

def test_init():
    def sequence_dummy(mol):
        return mol
    assert Standardizer()
    assert Standardizer(sequence_fun=sequence_dummy)
    assert Standardizer(sequence_fun=sequence_dummy, params=dict())

def test_sequence_minimal():
    # Violacein
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer().compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'
    assert MolToSmiles(ans) == 'OC1=NC(c2c[nH]c3ccc(O)cc23)=C/C1=C1\\C(O)=Nc2ccccc21'
    # L-Lactate
    mol = MolFromInchi('')

def test_sequence_rr_legacy():
    # Violacein
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_rr_legacy').compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'
    assert MolToSmiles(ans) == '[H]OC1=NC(C2=C([H])N([H])C3=C2C([H])=C(O[H])C([H])=C3[H])=C([H])/C1=C1\\C(O[H])=NC2=C([H])C([H])=C([H])C([H])=C21'

def test_sequence_tunable():
    # Check default arguments
    args, varargs, varkw, defaults, kwonlyargs, kwonlydefaults, annotations = inspect.getfullargspec(sequence_tunable)
    default_params = dict(zip(args[-len(defaults):], defaults))
    assert default_params == {
            'OP_REMOVE_ISOTOPE':True,
            'OP_NEUTRALISE_CHARGE': True,
            'OP_REMOVE_STEREO': False,
            'OP_COMMUTE_INCHI': False,
            'OP_KEEP_BIGGEST': True,
            'OP_ADD_HYDROGEN': True,
            'OP_KEKULIZE': True,
            'OP_NEUTRALISE_CHARGE_LATE': True
    }
    # Violacein, default parameter
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_tunable').compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'
    assert MolToSmiles(ans) == '[H]OC1=NC(C2=C([H])N([H])C3=C2C([H])=C(O[H])C([H])=C3[H])=C([H])/C1=C1\\C(O[H])=NC2=C([H])C([H])=C([H])C([H])=C21'
    # Violacein, strip stereo
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_REMOVE_STEREO': True}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)'
    assert MolToSmiles(ans) == '[H]OC1=C([H])C2=C(C([H])=C1[H])N([H])C([H])=C2C1=C([H])C(=C2C(=O)N([H])C3=C([H])C([H])=C([H])C([H])=C23)C(=O)N1[H]'
    # Violacien, implicit Hs
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_ADD_HYDROGEN': False}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'
    assert MolToSmiles(ans) == 'OC1=CC2=C(C=C1)NC=C2C1=C/C(=C2/C3=CC=CC=C3N=C2O)C(O)=N1'
    # Violacien, no kekulerization
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_KEKULIZE': False}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+'
    assert MolToSmiles(ans) == '[H]OC1=NC(c2c([H])n([H])c3c([H])c([H])c(O[H])c([H])c23)=C([H])/C1=C1\\C(O[H])=Nc2c([H])c([H])c([H])c([H])c21'
    # Violacien, strip stereo & implicit Hs & no kekulerization
    mol = MolFromInchi('InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)/b18-13+')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_REMOVE_STEREO': True, 'OP_ADD_HYDROGEN': False, 'OP_KEKULIZE': False}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C20H13N3O3/c24-10-5-6-15-12(7-10)14(9-21-15)17-8-13(19(25)23-17)18-11-3-1-2-4-16(11)22-20(18)26/h1-9,21,24H,(H,22,26)(H,23,25)'
    assert MolToSmiles(ans) == 'O=C1NC(c2c[nH]c3ccc(O)cc23)=CC1=C1C(=O)Nc2ccccc21'
    # Lactate, default parameter
    mol = MolFromSmiles('C[C@@H](C(=O)[O-])O')
    ans = Standardizer(sequence_fun='sequence_tunable').compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1'
    assert MolToSmiles(ans) == '[H]OC(=O)[C@@]([H])(O[H])C([H])([H])[H]'
    # L-lactate, implicit Hs
    mol = MolFromSmiles('C[C@@H](C(=O)[O-])O')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_ADD_HYDROGEN': False}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/t2-/m0/s1'
    assert MolToSmiles(ans) == 'C[C@H](O)C(=O)O'
    # L-lactate, no stereo
    mol = MolFromSmiles('C[C@@H](C(=O)[O-])O')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_REMOVE_STEREO': True}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
    assert MolToSmiles(ans) == '[H]OC(=O)C([H])(O[H])C([H])([H])[H]'
    # L-lactate, no charge neutralisation
    mol = MolFromSmiles('C[C@@H](C(=O)[O-])O')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_NEUTRALISE_CHARGE': False, 'OP_NEUTRALISE_CHARGE_LATE': False}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)/p-1/t2-/m0/s1'
    assert MolToSmiles(ans) == '[H]O[C@]([H])(C(=O)[O-])C([H])([H])[H]'
    # L-lactate, implicit Hs & no stereo
    mol = MolFromSmiles('C[C@@H](C(=O)[O-])O')
    ans = Standardizer(sequence_fun='sequence_tunable', params={'OP_ADD_HYDROGEN': False, 'OP_REMOVE_STEREO': True}).compute(mol)
    assert MolToInchi(ans) == 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
    assert MolToSmiles(ans) == 'CC(O)C(=O)O'
