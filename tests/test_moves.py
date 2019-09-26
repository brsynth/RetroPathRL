"""
Aim: test compound features
"""

# RP3 objects
from compound import Compound
from move import Move

class TestMove(object):
    """
    Testing moves - should be fast
    """
    def test_cloning(self):
        move = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id")
        cloned_move = move.clone()
        different_python_object = (id(move) != id(cloned_move))
        identical_move_object = move.eq_full_inchi_key(cloned_move)
        assert (different_python_object and identical_move_object)

    def test_equality_true(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]")
        move = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id",
                    product_list = [compound_1, compound_6],
                    set_number = 5)
        move_bis = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id",
                    product_list = [compound_6, compound_1])

        assert move.eq_full_inchi_key(move_bis)

    def test_equality_false(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        compound_2345 = Compound("[H][C](=[O])[C]([H])=[C]([H])[H]")
        move = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id",
                    product_list = [compound_1, compound_6])
        move_bis = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id",
                    product_list = [compound_6, compound_1, compound_2345])
        move_ter = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id_2",
                    product_list = [compound_6, compound_1])

        assert move != move_bis and move != move_ter and move_bis != move_ter

    def test_rave_update(self):
        compound_1 = Compound("[H+]")
        compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]")
        move = Move(rsmart = "rsmart",
                    rid = "rid",
                    compound_id= "compound_id",
                    product_list = [compound_1, compound_6])

        move.update(5, visit_number = 10)
        move.update(0.2, 10)
        assert move.RAVE_total_score == 52
        assert move.RAVE_visits == 20

    # def more_compelx_tests_wthi_compouns
