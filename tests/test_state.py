"""
Aim: test compound features
"""

# General utility packages
import random
import pickle

# RP3 specific packages
from compound import Compound, unpickle
from chemical_compounds_state import ChemicalCompoundState
from rewarding import Basic_Rollout_Reward
from organisms import detectable_cmpds_H, Test_organism_H
from organisms import detectable_cmpds_noH, Test_organism_noH
from rule_sets_examples import applicable_rules_mixed_dict

random.seed(42)

class TestState(object):
    """
    Testing compounds
    """
    def test_set_up_as_list_compounds(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState([compound])
        assert compound.eq_full_inchi_key(state.compound_list[0])

    def test_set_up_as_list_smiles(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState([csmile])
        print(state.compound_list[0])
        assert compound.eq_full_inchi_key(state.compound_list[0])

    def test_set_up_as_Compound(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState(compound)
        assert compound.eq_full_inchi_key(state.compound_list[0])

    def test_set_up_as_smiles(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState(csmile)
        assert compound.eq_full_inchi_key(state.compound_list[0])

    def test_cloning(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState(compound)
        cloned_state = state.clone()
        different_python_object = (id(state) != id(cloned_state))
        identical_chemical_object = (state == cloned_state)
        assert (different_python_object and identical_chemical_object)

    def test_saving_and_loading(self):
        compound_1 = Compound("[H+]")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]")
        chassis_1_10_909 = ChemicalCompoundState([compound_1, compound_10, compound_909], organism = Test_organism_H)
        chassis_1_10_909.save(folder_address = "tests/generated_jsons/pickles")
        chassis_load = unpickle("BOPG_BSAB_GPRL", type = "state", folder_address = "tests/data")
        assert chassis_1_10_909 == chassis_load

    def test_check_sanitized_state_not(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        unstandardised = "C(=O)CC"
        compound = Compound(csmile)
        compound_ns = Compound(unstandardised)
        state = ChemicalCompoundState([compound, compound_ns], main_layer = False)  # state is not sanitised
        assert not state._check_sanitized_state()

    def test_check_sanitized_state_not(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        unstandardised = "C(=O)CC"
        compound = Compound(csmile)
        compound_ns = Compound(unstandardised)
        state = ChemicalCompoundState([compound, compound_ns], main_layer = True)  # state is not sanitised
        assert state._check_sanitized_state()

    def test_check_sanitized_state_yes(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state = ChemicalCompoundState([compound, compound_2])  # state is not sanitised
        assert state._check_sanitized_state()

    def test_compound_in_state_yes(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        unstandardised = "C(=O)CC"
        compound = Compound(csmile)
        compound_ns = Compound(unstandardised)
        state = ChemicalCompoundState([compound], organism = Test_organism_H)  # state is not sanitised
        assert state.compound_in_state(compound_ns)

    # def test_compound_in_state_yes_bigger_state(self):
    #     csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
    #     csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
    #     unstandardised = "C(=O)CC"
    #     compound = Compound(csmile)
    #     compound_2 = Compound(csmile_2)
    #     compound_ns = Compound(unstandardised)
    #     state = ChemicalCompoundState([compound, compound_2])  # state is not sanitised
    #     assert state.compound_in_state(compound_ns)

    def test_compound_in_state_no(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state = ChemicalCompoundState([compound], organism = Test_organism_H)  # state is not sanitised
        assert not state.compound_in_state(compound_2)

    def test_add_compound_to_state_yes(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state = ChemicalCompoundState(compound, organism = Test_organism_H)
        state._add_compound_to_state(compound_2)
        assert len(state) == 2

    def test_add_compound_to_state_no(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile)
        state = ChemicalCompoundState(compound, organism = Test_organism_H)
        state._add_compound_to_state(compound_2)
        assert len(state) == 1

    def test_equal_two(slef):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state = ChemicalCompoundState([compound, compound_2], organism = Test_organism_H)
        state_2 = ChemicalCompoundState([compound, compound_2], organism = Test_organism_H)
        assert state == state_2

    def test_equal_two_from_addition(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state = ChemicalCompoundState(compound, organism = Test_organism_H)
        state_2 = ChemicalCompoundState(compound_2, organism = Test_organism_H)
        state._add_compound_to_state(compound_2)
        state_2._add_compound_to_state(compound)
        assert state == state_2

    def test_differing_compound_not_in_self(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]")
        state_1 = ChemicalCompoundState([compound_1, compound_6, compound_2345])
        state_2 = ChemicalCompoundState([compound_1, compound_6, compound_3459])
        assert state_1.not_in_self(state_2)[0] == compound_3459

    def test_differing_compound_not_in_other(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]")
        state_1 = ChemicalCompoundState([compound_1, compound_6, compound_2345])
        state_2 = ChemicalCompoundState([compound_1, compound_6, compound_3459])
        assert state_1.not_in_other(state_2)[0] == compound_2345

    # def test_GetMoves_archive(self):
    #     # will have to be modified once we have proper scores
    #     # random.seed(42)
    #     csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
    #     csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
    #     compound = Compound(csmile)
    #     compound_2 = Compound(csmile_2)
    #     state_2 = ChemicalCompoundState([compound, compound_2], organism = Test_organism_H)
    #     moves = state_2.GetMoves_archive(top_x = 2)
    #     assert list(moves.keys()) == ['1-MNXR94682_MNXM821', '1-MNXR117465_MNXM821']

    def test_GetMoves(self):
        # will have to be modified once we have proper scores
        # random.seed(42)
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        state_2 = ChemicalCompoundState([compound, compound_2], organism = Test_organism_H)
        moves = state_2.GetMoves(top_x = 2, chemical_score = False)
        for move in moves:
            print(move.rid)
        assert set([move.rid for move in moves]) == set(['MNXR103108_MNXM90191', 'MNXR95713_MNXM90191'])

    def test_GetMoves_state_and_compound(self):
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        state_test_10 = ChemicalCompoundState(compound_10, organism = Test_organism_H, available_rules = applicable_rules_mixed_dict)
        moves_compound = compound_10.obtain_applicable_transformation_with_move(available_rules = applicable_rules_mixed_dict)
        print("moves_compound: {}".format(moves_compound))
        del compound_10
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        state_test_10 = ChemicalCompoundState(compound_10, organism = Test_organism_H, available_rules = applicable_rules_mixed_dict)

        moves_state = state_test_10.GetMoves(top_x = 10)
        print("moves_state: {}".format(moves_state))
        identicial_moves = len(moves_state) == len(moves_compound)
        for move in moves_state:
            identicial_moves = identicial_moves and move.in_list(moves_compound, main_layer = False)
        assert identicial_moves
        # print(len(moves_compound))
        # # print("moves compounds {}".format(moves_compound))
        # print(len(moves_state))

    def test_ApplyMoves(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        compound_2 = Compound(csmile_2)
        known_result_state = ChemicalCompoundState([compound_2], organism = Test_organism_H)
        state_apply_move = ChemicalCompoundState([compound], organism = Test_organism_H)
        moves = state_apply_move.GetMoves(top_x = 2)
        # move = moves[list(moves.keys())[0]]
        for move in moves:
            if move.rid == "MNXR94682_MNXM821":
                i = moves.index(move)
        move = moves[i]
        state_apply_move.ApplyMove(move)
        assert state_apply_move == known_result_state

    def test_GetResults_from_Compounds_penalty(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState([compound], organism = Test_organism_H)
        assert state.GetResults_from_Compounds(rewarding = Basic_Rollout_Reward) == -1

    def test_GetResults_from_InChI_Key_penalty(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile)
        state = ChemicalCompoundState([compound], organism = Test_organism_H)
        assert state.GetResults_from_InChI_Keys(rewarding = Basic_Rollout_Reward) == -1

    def test_compound_in_chassis_state(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]")
        chassis_metabolites = ChemicalCompoundState([compound_1, compound_6, compound_3459], organism = Test_organism_H)
        assert chassis_metabolites.compound_in_state(compound_1)

    def test_GetResults_from_Compounds_reward(self):
        chassis_1_smile = "[H+]"
        chassis_2_smile = "[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]"
        chassis_1 = Compound(chassis_1_smile)
        chassis_2 = Compound(chassis_2_smile)
        state_chassis = ChemicalCompoundState([chassis_1, chassis_2], organism = Test_organism_H)
        assert state_chassis.GetResults_from_Compounds(rewarding = Basic_Rollout_Reward) == Basic_Rollout_Reward.full_state_reward

    def test_GetResults_from_InChI_Key_reward(self):
        chassis_1_smile = "[H+]"
        chassis_2_smile = "[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]"
        chassis_1 = Compound(chassis_1_smile)
        chassis_2 = Compound(chassis_2_smile)
        state_chassis = ChemicalCompoundState([chassis_1, chassis_2], organism = Test_organism_H)
        assert state_chassis.GetResults_from_InChI_Keys(rewarding = Basic_Rollout_Reward) == Basic_Rollout_Reward.full_state_reward

    def test_GetResults_from_Compounds_fraction(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        chassis_2_smile = "[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]"
        outside_compound = Compound(csmile)
        chassis_2 = Compound(chassis_2_smile)
        state_chassis = ChemicalCompoundState([chassis_2, outside_compound], organism = Test_organism_H)
        assert state_chassis.GetResults_from_Compounds(rewarding = Basic_Rollout_Reward) == 1/2

    def test_GetResults_from_InChI_Key_fraction(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        chassis_2_smile = "[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]"
        outside_compound = Compound(csmile)
        chassis_2 = Compound(chassis_2_smile)
        state_chassis = ChemicalCompoundState([chassis_2, outside_compound], organism = Test_organism_H)
        assert state_chassis.GetResults_from_InChI_Keys(rewarding = Basic_Rollout_Reward) == 1/2

    def test_GetResultsForBiosensors_present(self):
        chassis_1_smile = "[H+]"
        chassis_1 = Compound(chassis_1_smile)
        benzoate = Compound(InChI = "InChI=1S/C7H6O2/c8-7(9)6-4-2-1-3-5-6/h1-5H,(H,8,9)/p-1")
        benzoate_chassis = ChemicalCompoundState([benzoate, chassis_1], organism = detectable_cmpds_H)
        retrosynthesis = benzoate_chassis.GetResults_from_InChI_Keys(rewarding = Basic_Rollout_Reward)
        biosensor = benzoate_chassis.GetResultsForBiosensors(rewarding = Basic_Rollout_Reward)
        assert retrosynthesis ==0.5 and biosensor == 2

    def test_GetResultsForBiosensors_absent(self):
        chassis_1_smile = "[H+]"
        chassis_1 = Compound(chassis_1_smile)
        benzoate_chassis = ChemicalCompoundState([chassis_1], organism = detectable_cmpds_H)
        retrosynthesis = benzoate_chassis.GetResults_from_InChI_Keys(rewarding = Basic_Rollout_Reward)
        biosensor = benzoate_chassis.GetResultsForBiosensors(rewarding = Basic_Rollout_Reward)
        assert retrosynthesis ==-1 and biosensor == -1

    def test_remove_cmps(self):
        chassis_1_smile = "[H+]"
        chassis_2_smile = "[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]"
        chassis_1 = Compound(chassis_1_smile)
        chassis_2 = Compound(chassis_2_smile)
        state_chassis = ChemicalCompoundState([chassis_1, chassis_2], organism = Test_organism_H)
        present_in_state = state_chassis.compound_in_state(chassis_1)
        state_chassis.remove_cmpd_from_state(chassis_1)
        present_in_state_after = state_chassis.compound_in_state(chassis_1)
        assert present_in_state and not present_in_state_after

        # chassis_metabolites == ChemicalCompoundState([Compound("[H+]"),
        #     Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]"),
        #     Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]")])

    # Verify states are equal buy addign compoudns in 2 different ways
