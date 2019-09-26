"""
Aim: test compound features
"""

# General utility packages
import random
import pickle

# RP3 specific classes
from compound import Compound, unpickle
from chemical_compounds_state import ChemicalCompoundState
from MCTS_node import MCTS_node
from rewarding import Basic_Rollout_Reward
from organisms import Test_organism_H

random.seed(42)


class TestMCTSNode(object):

    def test_proper_root_initialisation(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        assert node._verif_initilisation()

    def test_terminal_state_false(self):
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound = Compound(csmile, name = "821")
        state = ChemicalCompoundState([compound], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        assert node.terminal == False

    def test_terminal_state_true_chassis(self):
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]", name = '3459')
        state = ChemicalCompoundState([compound_1, compound_6, compound_3459], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        assert node.terminal == True

    def test_terminal_state_true_unavailable(self):
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        state = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        assert node.terminal == True

    def test_terminal_state_true_unavailable_moves(self):
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        state = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        assert node.moves == []

    def test_proper_child_addition_821_rule_94682_moves(self):
        # random.seed(42)
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile, "821")
        compound_2 = Compound(csmile_2, "90191")
        known_result_state = ChemicalCompoundState([compound_2], organism = Test_organism_H)
        state_apply_move = ChemicalCompoundState(compound, organism = Test_organism_H)
        moves = state_apply_move.GetMoves(top_x = 2)
        # state_apply_move.ApplyMove(list(moves.keys())[0])
        assert set([move.rid for move in moves]) == set(['MNXR94682_MNXM821', 'MNXR117465_MNXM821'])

    def test_proper_child_addition_821_rule_94682_moves_state(self):
        # random.seed(42)
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile, "821")
        compound_2 = Compound(csmile_2, "90191")
        known_result_state = ChemicalCompoundState([compound_2], organism = Test_organism_H)
        state_apply_move = ChemicalCompoundState(compound, organism = Test_organism_H)
        moves = state_apply_move.GetMoves(top_x = 2)
        for move in moves:
            print(move.rid)
            if move.rid == "MNXR94682_MNXM821":
                i = moves.index(move)
        # print(moves)
        move_94682 = moves[i]
        state_apply_move.ApplyMove(move_94682)
        assert state_apply_move == known_result_state

    def test_proper_child_addition_821_rule_117465_moves(self):
        # random.seed(42)
        csmile = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_2 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile, "821")
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        known_result_state = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)
        state_apply_move = ChemicalCompoundState(compound, organism = Test_organism_H)
        moves = state_apply_move.GetMoves(top_x = 2)
        for move in moves:
            print(move.rid)
            if move.rid == "MNXR117465_MNXM821":
                i = moves.index(move)
        move_117465 = moves[i]
        state_apply_move.ApplyMove(move_117465)
        assert state_apply_move == known_result_state

    def test_expand_821(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        csmile_90191 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound = Compound(csmile_821, "821")
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        compound_90191 = Compound(csmile_90191, "90191")
        resulting_state_94682 = ChemicalCompoundState([compound_90191], organism = Test_organism_H)
        resulting_state_117465 = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)
        root_state = ChemicalCompoundState(compound, organism = Test_organism_H)
        root_node = MCTS_node(root_state, chemical_score = True)
        root_node.expand()
        for child in root_node.children:
            if child.move.rid == "MNXR94682_MNXM821":
                i_94682 = root_node.children.index(child)
            if child.move.rid == "MNXR117465_MNXM821":
                i_117465 = root_node.children.index(child)
        result_94682 = (root_node.children[i_94682].state == resulting_state_94682)
        result_117465 = (root_node.children[i_117465].state == resulting_state_117465)
        assert (result_94682 and result_117465)

    def test_expand_90191(self):
        csmile_90191 = "[H][O][C]([H])([H])[C]([H])([O][H])[C]([H])([H])[H]"
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = 909)
        # compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        compound_90191 = Compound(csmile_90191)
        resulting_state_95713 = ChemicalCompoundState([compound_1, compound_10, compound_909], organism = Test_organism_H)
        resulting_state_903108 = ChemicalCompoundState([compound_1, compound_6, compound_909], organism = Test_organism_H)
        root_state = ChemicalCompoundState(compound_90191, organism = Test_organism_H)
        root_node = MCTS_node(root_state)
        root_node.expand()

        result_95713 = (root_node.children[0].state == resulting_state_95713 or root_node.children[0].state == resulting_state_903108)
        result_903108 = (root_node.children[1].state == resulting_state_903108 or root_node.children[1].state == resulting_state_95713)
        assert (result_95713 and result_903108)

    def test_expand_1_6_2345_empty_children(self):
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        state = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        node.expand()
        assert node.children == []

    def test_expand_1_6_2345_same_node(self):
        compound_1 = Compound("[H+]", name = "1")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]", name = "2345")
        state = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)  # state is not sanitised
        node = MCTS_node(state)
        node.expand()
        state_ue = ChemicalCompoundState([compound_1, compound_6, compound_2345], organism = Test_organism_H)  # state is not sanitised
        node_ue = MCTS_node(state_ue)
        assert node == node_ue

    def test_node_equality_true(self):
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = 909)
        state_1 = ChemicalCompoundState([compound_6, compound_10, compound_909], organism = Test_organism_H)
        state_2 = ChemicalCompoundState([compound_6, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_2 = MCTS_node(state_2)
        assert (node_1 == node_2)

    def test_pickling_and_saving(self):
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = 909)
        state_1 = ChemicalCompoundState([compound_6, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_1.save(folder_address = "tests/generated_jsons/pickles")
        node_2 = unpickle("6_10_909-0-0", type = "node", folder_address = "tests/generated_jsons/pickles")
        assert (node_1 == node_2)

    def test_node_unexpanded(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = "909")
        state_1 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        assert (node_1.children == [])

    def test_node_expanded(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = "909")
        state_1 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_1.expand()
        assert (node_1.children != [])

    def test_node_equality_true_expanded(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = "909")
        state_1 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        state_2 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_2 = MCTS_node(state_2)
        node_2.expand()
        assert (node_1 == node_2)

    def test_node_equality_true_expanded(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = "909")
        state_1 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        state_2 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_2 = MCTS_node(state_2)
        node_2.expand()
        assert (node_1 == node_2)

    def test_node_equality_false_expanded_updated(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = "10")
        compound_909 = Compound("[H][O][C]([H])([C]([H])=[O])[C]([H])([H])[H]", name = "909")
        state_1 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        state_2 = ChemicalCompoundState([compound_821, compound_10, compound_909], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        node_2 = MCTS_node(state_2)
        node_2.update(10, solved = True)
        assert (node_1 != node_2)

    def test_node_rollout_possible(self):
        csmile_821 = "[H][C](=[O])[C]([H])([H])[C]([H])([H])[H]"
        compound_821 = Compound(csmile_821, "821")
        state_1 = ChemicalCompoundState([compound_821], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        moves = state_1.GetMoves(top_x = 2, chemical_score = False)
        rolled_state = node_1.rollout()
        assert (node_1.state != rolled_state)

    def test_node_rollout_terminal(self):
        compound_10 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", name = '6')
        state_1 = ChemicalCompoundState([compound_10, compound_6], organism = Test_organism_H)
        node_1 = MCTS_node(state_1)
        rolled_state = node_1.rollout()
        assert (rolled_state == node_1.state)
