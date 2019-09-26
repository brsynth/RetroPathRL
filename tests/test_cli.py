"""
Test FireBurner class
"""

import rdkit
from rdkit import Chem
import pytest
import multiprocessing


from utilities.reactor.cli import RuleBurner, RuleConversionError, ChemConversionError


# Data for tests
substate_inchi = 'InChI=1S/C3H6O3/c1-2(4)3(5)6/h2,4H,1H3,(H,5,6)'
reaction_smarts = '([#8&v2:1](-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11])-[#1&v1:12])>>([#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])(-[#8&v2]-[#1&v1])-[#8&v2:1]-[#6&v4:2](-[#6&v4:3](-[#8&v2:4]-[#1&v1:5])=[#8&v2:6])(-[#6&v4:7](-[#1&v1:8])(-[#1&v1:9])-[#1&v1:10])-[#1&v1:11].[#7&v3](=[#6&v4]1:[#7&v3]:[#6&v4](-[#8&v2]-[#1&v1]):[#6&v4]2:[#7&v3]:[#6&v4](-[#1&v1]):[#7&v3](-[#6&v4]3(-[#1&v1])-[#8&v2]-[#6&v4](-[#6&v4](-[#8&v2]-[#15&v5](=[#8&v2])(-[#8&v2]-[#1&v1])-[#8&v2]-[#15&v5](-[#8&v2]-[#1&v1:12])(=[#8&v2])-[#8&v2]-[#1&v1])(-[#1&v1])-[#1&v1])(-[#1&v1])-[#6&v4](-[#8&v2]-[#1&v1])(-[#1&v1])-[#6&v4]-3(-[#8&v2]-[#1&v1])-[#1&v1]):[#6&v4]:2:[#7&v3]:1-[#1&v1])-[#1&v1])'
tuple_product_inchikeys = ('CSZRNWHGZPKNKY-UHFFFAOYSA-N', 'QGWNDRXFNXRZMB-UHFFFAOYSA-N')
tuple_product_smiles = ('[H][O][C](=[O])[C]([H])([O][P](=[O])([O][H])[O][H])[C]([H])([H])[H]', '[H][N]=[c]1[n][c]([O][H])[c]2[n][c]([H])[n]([C]3([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][H])[C]([H])([O][H])[C]3([H])[O][H])[c]2[n]1[H]')
tuple_product_inchis = ('InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)', 'InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)')


def dummy_worker(**kwargs):
    import time
    time.sleep(1)


def test_init():
    # Empty is OK
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])  # Empty is OK
    rb.compute()


def test_run_with_timeout():
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])
    with pytest.raises(multiprocessing.context.TimeoutError):
        rb._run_with_timeout(dummy_worker, None, timeout=0)
    rb._run_with_timeout(dummy_worker, None, timeout=2)


def test_jsonify():
    rb = RuleBurner(rsmarts_list=[], inchi_list=[])
    assert rb._jsonify(rsmarts='', inchi='', rid='RID', cid='CID').replace('\n', '') == """{ "rule_id": "RID", "substrate_id": "CID", "fire_timed_out": null, "fire_exec_time": null}"""


def test_compute():
    # Wrong reaction depiction
    rb = RuleBurner(rsmarts_list=['DUMMY'], inchi_list=[])
    with pytest.raises(RuleConversionError):
        rb.compute()
    # Wrong chemical depiction
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=['DUMMY'])
    with pytest.raises(ChemConversionError):
        rb.compute()
    # Timeout should be logged
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=[substate_inchi], fire_timeout=0)
    rb.compute()
    assert ''.join(rb._json).find('"fire_timed_out": true')
    # OK
    rb = RuleBurner(rsmarts_list=[reaction_smarts], inchi_list=[substate_inchi])
    rb.compute()
    assert ''.join(rb._json).find('InChI=1S/C3H7O6P/c1-2(3(4)5)9-10(6,7)8/h2H,1H3,(H,4,5)(H2,6,7,8)')
    assert ''.join(rb._json).find('InChI=1S/C10H15N5O11P2/c11-10-13-7-4(8(18)14-10)12-2-15(7)9-6(17)5(16)3(25-9)1-24-28(22,23)26-27(19,20)21/h2-3,5-6,9,16-17H,1H2,(H,22,23)(H2,19,20,21)(H3,11,13,14,18)')
