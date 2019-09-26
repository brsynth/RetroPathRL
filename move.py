"""
Contains the move class that contains:
- compound it applies to
- rsmart
- rid
- set (becuase a signle rule cna generate multiple product sets)
- biological_score
- chemical_score
"""

# General utilities
import logging
import csv

class Move(object):
    """
    Basic move object. At the moment will have only attributes, no function per say.
    """

    logger = logging.getLogger(__name__)
    def __init__(self,
                 rsmart,
                 rid,
                 compound_id,
                 rsmiles = None,
                 set_number = 0,
                 chemical_score = 0,
                 chemical_substrate_score  = 0,
                 biological_score = 0,
                 product_list = [],
                 EC_number = ["EC: None"],
                 compound_index = 0, stoechiometry = {}):
        self.rsmart = rsmart
        if rsmiles is None:
            self.rsmiles = self.rsmart
        else:
            self.rsmiles = rsmiles
        self.rid = rid
        self.compound_id = compound_id
        self.set_number = set_number
        self.chemical_score = chemical_score
        self.chemical_substrate_score = chemical_substrate_score
        self.biological_score = biological_score
        self.EC_numbers = EC_number
        self.product_list = product_list
        self.name = "{}-{}-{}".format(self.compound_id, self.rid, str(self.set_number))
        self.synonyms = [self.rid]
        self.RAVE_visits = 0
        self.RAVE_total_score = 0
        self.RAVE_average_score = 0
        self.stoechiometry = stoechiometry

    def set_set_number(self, set_number):
        self.set_number = set_number
        self.name = "{}-{}-{}".format(self.compound_id, self.rid, str(self.set_number))

    def set_rsmart(self, rsmart):
        self.rsmart = rsmart

    def set_rsmiles(self, rsmiles):
        self.rsmiles = rsmiles

    def calculate_rsmiles(self, substrate):
        """
        Smiles of the actual transformation that is happening between the substrate and the products
        """
        sub_smiles = "{}".format(substrate.csmiles)
        prod_smiles = ".".join([prod.csmiles for prod in self.full_product_list()])
        self.rsmiles = "{}>>{}".format(sub_smiles, prod_smiles)

    def set_chemical_score(self, chemical_score):
        self.chemical_score = chemical_score

    def set_chemical_substrate_score(self, chemical_substrate_score):
        self.chemical_substrate_score = chemical_substrate_score

    def delete_intermediate_chemical_score(self):
        del self.original_substrates_list
        del self.original_products_list_list

    def set_intermediate_chemical_score(self, original_substrates_list, original_products_list_list):
        self.original_substrates_list = original_substrates_list
        self.original_products_list_list = original_products_list_list

    def set_id(self, id):
        self.id = id

    def set_EC_numbers(self, EC_numbers):
        self.EC_numbers = EC_numbers

    def set_biological_score(self, biological_score):
        self.biological_score = biological_score

    def set_product_list(self, product_list):
        self.product_list = product_list

    def set_stoechiometry(self, stoechiometry):
        self.stoechiometry = stoechiometry

    def __repr__(self):
        return(self.name)

    def print_all_attributes(self):
        text = "For move {}, attributes are: rid: {}, cid: {} \n".format(self.name, self.rid, self.compound_id)
        text_next = "set: {}, chem_score: {}, bio score: {} \n".format(self.set_number, self.chemical_score, self.biological_score)
        text_last = "product_list: {}, stoechiometry: {} \n".format(self.product_list, self.stoechiometry)
        text_appendix = "EC numbers are {}".format(self.EC_numbers)
        return(text + text_next + text_last + text_appendix)

    def full_product_list(self):
        full_list = []
        ordered_product_list = sorted(self.product_list, key = lambda item: self.stoechiometry[item.InChIKey])
        for product in ordered_product_list:
            for i in range(self.stoechiometry[product.InChIKey]):
                full_list.append(product)
        return(full_list)

    def _calculate_simles_from_move(self):
        sub = "{}".format()

    def clone(self):
        cloned_move = Move(rsmart = self.rsmart,
                        rid = self.rid,
                        compound_id = self.compound_id,
                        set_number = self.set_number,
                        chemical_score = self.chemical_score,
                        biological_score = self.biological_score,
                        product_list = self.product_list,
                        EC_number = self.EC_numbers,
                        stoechiometry = self.stoechiometry)
        try:
            cloned_move.set_intermediate_chemical_score(self.original_substrates_list, self.original_products_list_list)
        except AttributeError:
            pass
        return(cloned_move)

    def add_synonym(self, move):
        """
        Adds a synonym to this move.
        (When another move was deemed equal to current move (self))
        """
        if move.rid not in self.synonyms:
            self.synonyms.append(move.rid)
            for EC in move.EC_numbers:
                if EC not in self.EC_numbers:
                    self.EC_numbers.append(EC)
            if self.biological_score * self.chemical_score < move.biological_score * move.chemical_score:
                self.biological_score = move.biological_score
                self.chemical_score = move.chemical_score
                self.stoechiometry = move.stoechiometry

    def eq_full_inchi_key(self, other):
        """
        Tow moves are identical if they
        - apply to the same compound
        - generate the same products
        """
        compound_eq = (self.compound_id == other.compound_id)
        products_eq = len(self.product_list) == len(other.product_list)
        for product in self.product_list:
            products_eq = products_eq and (product.in_list(other.product_list, main_layer = False))
        return(compound_eq and products_eq)

    def eq_main_layer(self, other):
        """
        Tow moves are identical if they
        - apply to the same compound
        - generate the same products
        """
        compound_eq = (self.compound_id == other.compound_id)
        products_eq = len(self.product_list) == len(other.product_list)
        for product in self.product_list:
            products_eq = products_eq and (product.in_list(other.product_list, main_layer = True))
        return(compound_eq and products_eq)

    def in_list(self, list_moves, main_layer = False):
        in_list = False
        for move_in_list in list_moves:
            if main_layer:
                equality = self.eq_main_layer(move_in_list)
                if equality:
                    in_list = True
                    move_in_list.add_synonym(self)
                    break
            else:
                equality = self.eq_full_inchi_key(move_in_list)
                if equality:
                    in_list = True
                    move_in_list.add_synonym(self)
                    break
        return(in_list)

    def update(self, result, visit_number = 1):
        """
        Values are used only for RAVE implementation.
        """
        self.RAVE_visits = self.RAVE_visits + visit_number
        self.RAVE_total_score = self.RAVE_total_score + result * visit_number
        self.RAVE_average_score = self.RAVE_total_score/self.RAVE_visits
