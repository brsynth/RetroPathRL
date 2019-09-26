"""
Aim: test compound features
"""

# General utility packages
import random
import pickle

# RP3  specific objects
from compound import Compound, unpickle
from chemical_compounds_state import ChemicalCompoundState
from representation import Test_representation, Test_to_file
from organisms import detectable_cmpds_H, Test_organism_H
from organisms import detectable_cmpds_noH
from rewarding import Basic_Rollout_Reward
from MCTS_node import MCTS_node
from UCT_policies import Biochemical_UCT_1, Nature_UCT, Classical_UCT_RAVE, Classical_UCT_with_bias, Classical_UCT
from rule_sets_examples import applicable_rules_mixed_dict, applicable_rules_10_dict
from Tree import Tree
from rule_sets_similarity import get_rules_and_score, full_rules_forward_H, full_rules_retro_H, full_rules_forward_no_H, full_rules_retro_no_H



random.seed(42)


class TestTree(object):
    def test_equality_statement_not_expanded(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 100)
        assert test_Tree ==  test_Tree_bis

    def test_equality_statement_expanded(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 100)
        test_Tree.run_search()
        test_Tree_bis.run_search()
        assert test_Tree ==  test_Tree_bis

    def test_equality_statement_expanded_differnet_iter(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 1000)
        test_Tree.run_search()
        test_Tree_bis.run_search()
        assert test_Tree != test_Tree_bis

    def test_equality_statement_expanded_false(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 100)
        test_Tree.run_search()
        assert test_Tree != test_Tree_bis

    def test_equality_statement_expanded_states(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100, available_rules = applicable_rules_mixed_dict)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 500, available_rules = applicable_rules_mixed_dict)
        test_Tree.run_search()
        test_Tree_bis.run_search()
        different_trees = test_Tree != test_Tree_bis
        same_states = test_Tree.equality_visited_states(test_Tree_bis)
        assert different_trees and same_states

    def test_equality_statement_expanded_states_other_policies(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 100)
        test_Tree_bis = Tree(root_state = state_bis, itermax = 1000, UCT_policy = "Nature_UCT")
        test_Tree.run_search()
        test_Tree_bis.run_search()
        different_trees = test_Tree != test_Tree_bis
        same_states = test_Tree.equality_visited_states(test_Tree_bis)
        assert different_trees and same_states

    def test_pickling_unpickling(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 10000,  parallel = False,
                         Rollout_policy = "Rollout_policy_first",
                         UCT_policy = "Biochemical_UCT_1")
        test_Tree.run_search()
        test_Tree.save("test", folder_address = "tests/generated_jsons/pickles")
        loaded_tree = unpickle(file_name = 'test', type = 'tree', folder_address = "tests/generated_jsons/pickles")
        assert test_Tree == loaded_tree

    def test_pickling_unpickling_differ(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised
        state_bis = ChemicalCompoundState([compound], organism = Test_organism_H, representation = Test_representation)  # state is not sanitised

        test_Tree = Tree(root_state = state, itermax = 10000,  parallel = False,
                         Rollout_policy = "Rollout_policy_first",
                         UCT_policy = "Biochemical_UCT_1")
        test_Tree.run_search()
        test_Tree.save("test", folder_address = "tests/generated_jsons/pickles")
        test_Tree.run_search()
        loaded_tree = unpickle(file_name = 'test', type = 'tree', folder_address = "tests/generated_jsons/pickles")
        assert test_Tree != loaded_tree

    def test_biosensor(self):
        organism = detectable_cmpds_H
        inchi = "InChI=1S/C6H11NO2/c8-6(9)5-3-1-2-4-7-5/h5,7H,1-4H2,(H,8,9)"
        compound = Compound(InChI = inchi, name = "pipecolate")
        present_in_state_detectable = organism.compound_in_state(compound)
        if present_in_state_detectable:
            logging.warning("Removed compound from the detectable set to force enzymatic detection")
            organism.remove_cmpd_from_state(compound)
        rules, biological_scoring = get_rules_and_score(full_rules_forward_H = full_rules_forward_H,
                                                    full_rules_retro_H = full_rules_retro_H,
                                                    full_rules_forward_no_H = full_rules_forward_no_H,
                                                    full_rules_retro_no_H = full_rules_retro_no_H,
                                                    add_Hs = True,
                                                    retro = False,
                                                    diameters = [10, 12, 14, 16],
                                                    small = False,
                                                    c_name = None,
                                                    filtering_EC = ["1.5.3.7", "1.5.3"])
        state = ChemicalCompoundState([compound])  # state is not sanitised
        test_Tree = Tree(root_state = state, itermax = 1000,  parallel = False,
                         Rollout_policy = "Rollout_policy_first",
                         UCT_policy = "Biochemical_UCT_1", available_rules = rules, organism  = organism,
                         biological_scorer = biological_scoring,
                         folder_to_save = "tests/generated_jsons")
        test_Tree.run_search()
        loaded_tree = unpickle(file_name = 'pipecolate_test', type = 'tree', folder_address = "tests/data")
        same_states = test_Tree.equality_visited_states(loaded_tree)
        assert same_states
