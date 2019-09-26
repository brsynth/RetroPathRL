"""
Aim: test compound features
"""

# General utility packages
import pickle

# Chemistry packages
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
# from reactor.Core import RuleBurnerCore
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError
from utilities.reactor.cli import worker_match, worker_fire, kill

# Class to test
from compound import Compound, unpickle
# Scoring packages - RP3 objects
from biological_scoring import RandomBiologicalScorer
from chemical_scoring import RandomChemicalScorer, SubstrateChemicalScorer, SubandprodChemicalScorer

# Chemcial rule handlers - RP3 objects
from move import Move
from rule_sets_examples import applicable_rules_10_dict
from rule_sets_similarity import full_rules_retro_H, applicable_rules_mixed_dict_sim, applicable_rules_10_dict_sim

# General script configuration
from config import *

class TestCompound(object):
    """
    Testing compounds
    """
    def test_set_up(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        assert csmile == Chem.MolToSmiles(compound.rdmol,isomericSmiles=True, allHsExplicit=True)

    def test_simple_pickle(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound.save("test", folder_address = "tests/generated_jsons/pickles")
        loaded = unpickle("test", folder_address = "tests/generated_jsons/pickles")
        assert compound.eq_full_inchi_key(loaded)

    def test_standardisation(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        unstandardised = "C(=O)CC"
        ns_compound = Compound(unstandardised)
        assert csmile == Chem.MolToSmiles(ns_compound.rdmol,isomericSmiles=True, allHsExplicit=True)

    def test_standardisation_csmiles_attribute(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        unstandardised = "C(=O)CC"
        ns_compound = Compound(unstandardised)
        s_compound = Compound(csmile)
        assert ns_compound.csmiles == s_compound.csmiles

    def test_unstandardised_pickle(self):
        csmile = "C(=O)CC"
        compound = Compound(csmile)
        compound.save("test", folder_address = "tests/generated_jsons/pickles")
        loaded = unpickle("test", folder_address = "tests/generated_jsons/pickles")
        assert compound.eq_full_inchi_key(loaded)

    def test_InChI_keys(self):
        InChI_key = "BQJCRHHNABKAKU-UHFFFAOYSA-N"
        # Remark: we are not generating standard InChI keys or InChIs.
        # It's okay as long as we're consistent
        canonical_smiles = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
        morphine = Compound(canonical_smiles)
        assert morphine.InChIKey == InChI_key

    def test_name_none_given(self):
        InChI_key = "BQJCRHHNABKAKU"
        # Remark: we are not generating standard InChI keys or InChIs.
        # It's okay as long as we're consistent
        canonical_smiles = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
        morphine = Compound(canonical_smiles)
        assert morphine.name.split("-")[0] == InChI_key

    def test_name_morphine_given(self):
        InChI_key = "BQJCRHHNABKAKU-UHFFFAOYSA-N"
        # Remark: we are not generating standard InChI keys or InChIs.
        # It's okay as long as we're consistent
        canonical_smiles = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
        morphine = Compound(canonical_smiles, name = "morphine")
        assert morphine.name == 'morphine'

    def test_naming(self):
        InChI_key = "BQJCRHHNABKAKU-UHFFFAOYSA-N"
        # Remark: we are not generating standard InChI keys or InChIs.
        # It's okay as long as we're consistent
        canonical_smiles = "CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O"
        morphine = Compound(canonical_smiles)
        morphine.naming("morphine")
        assert morphine.name == 'morphine'

    def test_kekularisation(self):
        pass

    def test_stereo(self):
        pass

    def test_equality_main_layer(self):
        # TODO: define a proper test for this one (not passing full layer but passing main)
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_bis = Compound(csmile)
        assert compound.eq_main_layer(compound_bis)

    def test_equality_full_layer(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_bis = Compound(csmile)
        assert compound.eq_full_inchi_key(compound_bis)

    def test_cloning(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        cloned_compound = compound.clone()
        different_python_object = (id(compound) != id(cloned_compound))
        identical_chemical_object = (compound.eq_full_inchi_key(cloned_compound))
        assert (different_python_object and identical_chemical_object)

    def test_rdmol_import(self):
        ethylene_glycol = Compound(InChI = "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2")
        mol_file = Chem.inchi.MolFromInchi("InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2", sanitize=False)
        mol_file = standardize_chemical(mol_file)
        ethylene_glycol_bis = Compound(rdkit_obj = mol_file)
        assert(ethylene_glycol.eq_full_inchi_key(ethylene_glycol_bis))

    # def test_rule_matcher(self):  # aslo add more for stereo if we decide to do taht
    #     one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
    #     one_set_one_product_rsmart = "([#6&v4:1](=[#8&v2:2])(-[#6&v4:3](-[#6&v4:4](-[#1&v1:5])(-[#1&v1:6])-[#1&v1:7])(-[#1&v1:8])-[#1&v1:9])-[#1&v1:10])>>([#6&v4:1](-[#8&v2:2]-[#1&v1:9])(-[#6&v4:3](-[#6&v4:4](-[#1&v1:5])(-[#1&v1:6])-[#1&v1:7])(-[#8&v2]-[#1&v1])-[#1&v1:8])(-[#1&v1:10])-[#1&v1])"
    #     compound = Compound(one_set_one_product_compound)
    #     assert compound._match_transformation(rsmarts = one_set_one_product_rsmart)

    # def test_apply_transformation_one_set_one_product(self):
    #     one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
    #     one_set_one_product_rsmart = "([#6&v4:1](=[#8&v2:2])(-[#6&v4:3](-[#6&v4:4](-[#1&v1:5])(-[#1&v1:6])-[#1&v1:7])(-[#1&v1:8])-[#1&v1:9])-[#1&v1:10])>>([#6&v4:1](-[#8&v2:2]-[#1&v1:9])(-[#6&v4:3](-[#6&v4:4](-[#1&v1:5])(-[#1&v1:6])-[#1&v1:7])(-[#8&v2]-[#1&v1])-[#1&v1:8])(-[#1&v1:10])-[#1&v1])"
    #     one_set_one_product_smile = ['[H]OC([H])([H])C([H])(O[H])C([H])([H])[H]']
    #     compound = Compound(one_set_one_product_compound)
    #     list_list_smiles = compound.apply_transformation_move(rsmarts = one_set_one_product_rsmart, available_rules = applicable_rules_10_dict)
    #     for result_set in list_list_smiles:
    #         for compound in result_set:
    #             assert compound in one_set_one_product_smile

    # def test_apply_transformation_2_results(self):
    #     one_set_two_products_compound = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
    #     one_set_two_products_rsmart = "([#6&v4:1](-[#8&v2:2]-[#1&v1:3])(-[#6&v4:4](-[#6&v4:5](-[#1&v1:6])(-[#1&v1:7])-[#1&v1:8])(-[#8&v2:9]-[#1&v1:10])-[#1&v1:11])(-[#1&v1:12])-[#1&v1:13])>>([#6&v4:1](=[#8&v2:2])(-[#6&v4:4](-[#6&v4:5](-[#1&v1:6])(-[#1&v1:7])-[#1&v1:8])(-[#1&v1:11])-[#1&v1:3])-[#1&v1:12].[#8&v2:9](-[#1&v1:10])-[#1&v1:13])"
    #     one_set_two_products_smile = ["[H]O[H]", "[H]C(=O)C([H])([H])C([H])([H])[H]"]
    #     compound = Compound(one_set_two_products_compound)
    #     list_list_smiles = compound.apply_transformation_archive(rsmarts = one_set_two_products_rsmart, available_rules = applicable_rules_10_dict)
    #     for result_set in list_list_smiles:
    #         for compound in result_set:
    #             assert compound in one_set_two_products_smile
        # One with only 1 Product

    def test_apply_transformation_sets(self):

        set_1 = [Compound('[H]OC([H])([H])C([H])(O[H])C([H])([H])[H]').name]
        set_2 = [Compound('[H]C(=O)C([H])=C([H])[H]').name, Compound('[H+]').name, Compound('[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]').name]
        set_test_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(set_test_compound)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        found_set_1 = False
        found_set_2 = False
        all_found = True
        for move in moves:
            products = compound.apply_transformation_with_move(move)
            if set([product.name for product in products]) == set(set_1):
                found_set_1 = True
            elif set([product.name for product in products]) == set(set_2):
                found_set_2 = True
            else:
                all_found = False

        assert (found_set_1 and found_set_2 and len(moves) == 2)

    def test_apply_transformation_sets_parallel(self):

        set_1 = [Compound('[H]OC([H])([H])C([H])(O[H])C([H])([H])[H]').name]
        set_2 = [Compound('[H]C(=O)C([H])=C([H])[H]').name, Compound('[H+]').name, Compound('[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]').name]
        set_test_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(set_test_compound, single_worker = True)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        found_set_1 = False
        found_set_2 = False
        all_found = True
        for move in moves:
            products = compound.apply_transformation_with_move(move)
            if set([product.name for product in products]) == set(set_1):
                found_set_1 = True
            elif set([product.name for product in products]) == set(set_2):
                found_set_2 = True
            else:
                all_found = False

        assert (found_set_1 and found_set_2 and len(moves) == 2)

    def test_apply_transformation_sets_different_products(self):

        # set_1 = ['[H]OC([H])([H])C([H])(O[H])C([H])([H])[H]']
        # set_2 = ['[H]C(=O)C([H])=C([H])[H]', '[H+]', '[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]']
        set_test_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(set_test_compound)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        move_1 = moves[0]
        product_set_1 = compound.apply_transformation_with_move(move_1)
        move_2 = moves[1]
        product_set_2 = compound.apply_transformation_with_move(move_2)
        assert set([product.name for product in product_set_1]) != set([product.name for product in product_set_2])

    def test_apply_transformation_sets_different_products_parallel(self):

        # set_1 = ['[H]OC([H])([H])C([H])(O[H])C([H])([H])[H]']
        # set_2 = ['[H]C(=O)C([H])=C([H])[H]', '[H+]', '[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]']
        set_test_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(set_test_compound, single_worker = True)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        move_1 = moves[0]
        product_set_1 = compound.apply_transformation_with_move(move_1)
        move_2 = moves[1]
        product_set_2 = compound.apply_transformation_with_move(move_2)
        assert set([product.name for product in product_set_1]) != set([product.name for product in product_set_2])

    # def test_obtain_applicable_transformation(self):
    #     one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
    #     one_set_one_product_matches = ['MNXR94682_MNXM821', 'MNXR117465_MNXM821']
    #     compound = Compound(one_set_one_product_compound)
    #     match_IDS, match_dict = compound.obtain_applicable_transformation_archive(available_rules = applicable_rules_10_dict)
    #     print(match_IDS)
    #     # issue with rule and reaction IDs
    #     for match in match_IDS:
    #         assert match in one_set_one_product_matches

    def test_obtain_applicable_transformation_with_moves(self):
        one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        one_set_one_product_matches = ['MNXR94682_MNXM821', 'MNXR117465_MNXM821']
        compound = Compound(one_set_one_product_compound)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        # issue with rule and reaction IDs
        for move in moves:
            assert move.rid in one_set_one_product_matches

    # def test_obtain_applicable_transformation_none(self):
    #     unmatched_substrate = "[H][C](=[O])[C]([H])=[C]([H])[H]"
    #     compound = Compound(unmatched_substrate)
    #     match_IDS, match_dict = compound.obtain_applicable_transformation_archive(available_rules = applicable_rules_10_dict)
    #     assert match_IDS == []

    def test_obtain_applicable_transformation_with_move_none(self):
        unmatched_substrate = "[H][C](=[O])[C]([H])=[C]([H])[H]"
        compound = Compound(unmatched_substrate)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict)
        assert moves == []

    def test_similarity_substrate_only_equal_1(self):
        one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(one_set_one_product_compound)
        move_id = "MNXR94682_MNXM821"
        original_substrates_list = applicable_rules_10_dict_sim[move_id]["substrate_ECFP"]
        scorer = SubstrateChemicalScorer
        score = scorer.calculate(compound, original_substrates_list = original_substrates_list,
        original_products_list_list = None)
        assert score == 1

    def test_similarity_substrate_only_different(self):
        one_set_one_product_compound = "[H][O][C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(one_set_one_product_compound)
        move_id = "MNXR94682_MNXM821"
        original_substrates_list = applicable_rules_10_dict_sim[move_id]["substrate_ECFP"]
        scorer = SubstrateChemicalScorer
        score = scorer.calculate(compound, original_substrates_list = original_substrates_list,
        original_products_list_list = None)

        assert round(score, 2) == 0.21

    def test_similarity_substrate_and_product_equal_1(self):
        one_set_one_product_compound = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(one_set_one_product_compound)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_10_dict_sim)
        move = moves[1]
        original_substrates_list = applicable_rules_10_dict_sim[move.rid]["substrate_ECFP"]
        original_product_list = applicable_rules_10_dict_sim[move.rid]["products_ECFP"]
        scorer = SubandprodChemicalScorer
        score = scorer.calculate(compound,
                            products = move.product_list,
                            original_substrates_list = original_substrates_list,
                            original_products_list_list = original_product_list)
        assert score == 1

    def test_similarity_substrate_and_product_different_1(self):
        one_set_one_product_compound = "[H]N=C(O[H])C1=C([H])N(C2([H])OC([H])(C([H])([H])OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C3([H])OC([H])(n4c([H])nc5c(N([H])[H])nc([H])nc54)C([H])(OP(=O)(O[H])O[H])C3([H])O[H])C([H])(O[H])C2([H])O[H])C([H])=C([H])C1([H])[H]"

        compound = Compound(one_set_one_product_compound)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = applicable_rules_mixed_dict_sim)
        move = moves[0]
        original_substrates_list = applicable_rules_mixed_dict_sim[move.rid]["substrate_ECFP"]
        original_product_list = applicable_rules_mixed_dict_sim[move.rid]["products_ECFP"]
        scorer = SubandprodChemicalScorer
        score = scorer.calculate(compound,
                            products = move.product_list,
                            original_substrates_list = original_substrates_list,
                            original_products_list_list = original_product_list)
        assert round(score, 2) == 0.24

    # def test_stoechiometry(self):
    #     one_set_one_product_compound = "InChI=1S/C40H56O2/c1-31(2)17-13-20-34(5)23-15-25-35(6)24-14-21-32(3)18-11-12-19-33(4)22-16-26-36(7)27-28-40-38(8,9)29-37(41)30-39(40,10)42-40/h11-12,14-19,21-28,37,41H,13,20,29-30H2,1-10H3/b12-11+,21-14+,22-16+,25-15+,28-27+,32-18+,33-19+,34-23+,35-24+,36-26+"
    #     compound = Compound(InChI = one_set_one_product_compound, max_moves = 4)
    #     moves = compound.obtain_applicable_transformation_with_move(available_rules = full_rules_12)
    #     len_1 = len(moves[0].full_product_list()) == 3 and len(moves[0].product_list) == 2
    #     len_2 = len(moves[1].full_product_list()) == 1 and len(moves[1].product_list) == 1
    #     assert len_1 and len_2

    def test_stoechiometry_3(self):
        ongoing_cmp_inchi = "InChI=1S/C15H12O4/c16-11-8-13(18)15(14(19)9-11)12(17)7-6-10-4-2-1-3-5-10/h1-9,16,18-19H"
        compound = Compound(InChI = ongoing_cmp_inchi, max_moves = 5)
        moves = compound.obtain_applicable_transformation_with_move(available_rules = full_rules_retro_H, chemical_scorer = "SubandprodChemicalScorer")
        move_1_stoechio = moves[4].stoechiometry == {'LTYOQGRJFJAKNA-UHFFFAOYSA-N': 3, 'JVNVHNHITFVWIX-UHFFFAOYSA-N': 1}
        move_2_stoechio = moves[0].stoechiometry == {'JFRPZSBKZZVAHT-UHFFFAOYSA-N': 1}
        assert move_1_stoechio and move_2_stoechio

    def test_rule_with_timeout(self):
        ongoing_cmp_inchi = "InChI=1S/C53H78O4/c1-41(2)19-11-20-42(3)21-12-22-43(4)23-13-24-44(5)25-14-26-45(6)27-15-28-46(7)29-16-30-47(8)31-17-32-48(9)33-18-34-49(10)39-40-57-53(56)51-37-35-50(36-38-51)52(54)55/h19,21,23,25,27,29,31,33,35-39H,11-18,20,22,24,26,28,30,32,34,40H2,1-10H3,(H,54,55)"
        compound = Compound(InChI = ongoing_cmp_inchi, max_moves = 10)
        move_rid_list = ["RR-02-dbef5e1eacd20afc-06-F",
                        "RR-02-98f27ab4082414aa-02-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F",
                        "RR-02-fa79f7a2d3a06ca5-04-F"]
        moves = compound.obtain_applicable_transformation_with_move(available_rules = full_rules_retro_H, chemical_scorer = "SubandprodChemicalScorer")
        all_found = True
        for move in moves:
            found = move.rid in move_rid_list
            all_found = all_found and found
        assert all_found
