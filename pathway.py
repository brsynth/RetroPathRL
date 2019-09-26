"""
Contains the pathway objects for visualisation and export
"""

# General utilities
import logging
import csv
import copy
import json
import pickle
# RP3 specific objects
from compound import Compound
from move import Move
from chemical_compounds_state import ChemicalCompoundState
from organisms import Test_organism_H


class Pathway(object):
    """
    Pathway object.
    Has methods for quick visualisation as well as export to json (for visualisation and treatment)
    Also has cloning and compound addition
    """
    logger = logging.getLogger(__name__)

    def __init__(self, first_iteration = -1, target = None, compounds = [], moves = [],
                 file_to_save = "temporary_pathway_json", main_layer = True,
                 organism = Test_organism_H, edges = [], nodes_compounds = [], nodes_transformations = []):
        """
        Initialising a pathway object.
        A compound has an ID and a dict with chemical structures
        A reaction links 2 compounds and has a smart, scores etc
        self.compounds is a dictionnary of ID: chemical_struct_of_compound
        Remarks:
        - a pathway can only be defined for a fully solved Node (ie: in the Tree, not in rollout)
        - it needs to verify at each step what products are formed
        as those could have been deleted in the tree search (already in state)
        """
        self.first_iteration = first_iteration
        self.target = target
        self.organism = organism
        self.main_layer = main_layer
        self.compounds = compounds
        self.moves = moves
        self.file_to_save = file_to_save
        self.nodes_compounds = nodes_compounds
        self.nodes_transformations = nodes_transformations
        self.edges = edges
        self.pathway_as_dict = None

    def __eq__(self, other):
        """
        Two pathways are identical if their compounds and moves are identical
        """
        node_compounds_equal = len(self.nodes_compounds) == len(other.nodes_compounds)
        node_trasnfo_equal = len(self.nodes_transformations) == len(other.nodes_transformations)
        node_edges_equal = len(self.edges) == len(other.edges)
        compounds_equal = len(self.compounds) == len(other.compounds)
        if compounds_equal:
            for compound in self.compounds:
                in_other = compound.in_list(other.compounds, main_layer = True)
                if not in_other:
                    compounds_equal = False
                    break
        moves_equal = len(self.moves) == len(other.moves)
        if moves_equal:
            for move in self.moves:
                in_other = move.in_list(other.moves, main_layer = True)
                if not in_other:
                    moves_equal = False
                    break
        equality = compounds_equal and moves_equal and node_compounds_equal and node_trasnfo_equal and node_edges_equal
        return (equality)

    def __repr__(self):
        """
        Print list of compoudns and list of moves
        """
        rep = 'Compound \n'
        for compound in self.compounds:
            rep = rep + str(compound) + "\n"
        rep = rep + 'Edges \n'
        for edge in self.edges:
            rep = rep + edge["data"]["id"] + "\n"
        return(rep)

    def all_attributes_with_nodes(self):
        """
        Print list of compounds and list of moves
        """
        rep = 'Compound \n'
        for compound in self.compounds:
            rep = rep + str(compound) + "\n"
        rep = rep + 'Edges \n'
        for edge in self.edges:
            rep = rep + edge["data"]["id"] + "\n"
        for node_cp in self.nodes_compounds:
            rep = rep + node_cp["data"]["id"] + "\n"
        for node_tf in self.nodes_transformations:
            rep = rep + node_tf["data"]["id"] + "\n"
        return(rep)

    def set_file_to_save(self, file_to_save):
        self.file_to_save = file_to_save

    def set_main_layer(self, main_layer):
        self.main_layer = main_layer

    def set_first_iteration(self, first_iteration):
        self.first_iteration = first_iteration

    def clone(self):
        """ Cloning """
        duplicated_pathway = Pathway(first_iteration = self.first_iteration,
                                     organism = self.organism ,
                                     main_layer = self.main_layer,
                                     target = self.target,
                                     compounds = [cmp.clone() for cmp in self.compounds],
                                     moves = [mv.clone() for mv in self.moves],
                                     edges = copy.deepcopy(self.edges),
                                     nodes_compounds = copy.deepcopy(self.nodes_compounds),
                                     nodes_transformations = copy.deepcopy(self.nodes_transformations)
                                     )
        return(duplicated_pathway)

    def save(self, file_name = None, folder_address = "pickled_data"):
        if file_name is None:
            base_name = self.file_to_save
        file_saving = open('{}/pathway_{}.pkl'.format(folder_address, file_name), 'wb')
        pickle.dump(self, file_saving)

    def add_compound(self, compound, in_sink = None, is_source = 0):
        """
        Adding a compound object to the pathway.
        """
        if is_source:
            self.target = compound
        if not compound.in_list(self.compounds, main_layer = self.main_layer):
            self.compounds.append(compound)
            if in_sink is None:
                if self.organism.compound_in_state(compound):
                    in_sink = 1
                else:
                    in_sink = 0
            data_dict = {
                'SMILES': compound.csmiles,
                'inSink':in_sink,
                'isSource': is_source,
                'InChI': compound.InChI,
                'Names': compound.synonyms_names,  # If I want synonyms, keep them
                'id': compound.InChIKey,
                'type': 'compound',
                'Rule ID': None,
                'EC number': None,
                'Reaction SMILES': None,
                'Diameter': None,
                'Score': None,
                'Iteration': None
                }
            self.nodes_compounds.append({"data": data_dict})
        else:
            self.logger.warning("Compound {} is already in compounds".format(compound))

    def clean_up(self, move, depth):
        str = "{}-{}-{}-{}".format(move.compound_id, move.rid, move.set_number, depth)
        return(str)

    def add_reaction(self, move, depth = 1):
        """
        Adding a reaction to the pathway.
        """
        if not move.in_list(self.moves):
            self.moves.append(move)
            move_compound_id_present = False
            for cp in self.compounds:
                for sym in cp.synonyms_names:
                    if sym == move.compound_id:
                        move_compound_id_present = True
                        move_compound_ID = cp.InChIKey
                        break
            if not move_compound_id_present:
                self.logger.warning("Trying to add move {} when compound {} is not in the pathway".format(move, move.compound_id))

            for product in move.product_list:
                if not product.in_list(self.compounds):
                    # Adding the products of the pathway
                    self.add_compound(product, in_sink = None, is_source = 0)

            cleaned_up_moved = self.clean_up(move, depth)
            try:
                diameter = int(move.rid.split("-")[3])
            except:
                diameter = 42
            data_dict = {
                "SMILES": None,
                "inSink": None,
                "isSource": None,
                "InChI": None,
                "Names": None,
                "id": cleaned_up_moved,
                "type": "reaction",
                "Rule ID": move.synonyms,
                "EC number": move.EC_numbers,
                "Reaction SMILES": move.rsmiles,
                "Diameter": diameter,
                "Score": move.biological_score,
                "ChemicalScore": move.chemical_score,
                "Iteration": depth,
                "Stoechiometry": move.stoechiometry
            }
            self.nodes_transformations.append({"data": data_dict})
            # Adding all the edges:
            # from compound to reaction (move as target, compound as source)
            # From reactions to compound (move as source, product as target)
            data_dict = {
                "target" : cleaned_up_moved,
                "source" : move_compound_ID,
                "id" : "{}_=>_{}".format(cleaned_up_moved, move.compound_id)
            }
            self.edges.append({"data": data_dict})
            for product in move.product_list:
                data_dict = {
                    "target" : product.name,
                    "source" : cleaned_up_moved,
                    "id" : "{}_=>_{}".format(product.name, cleaned_up_moved)
                }
                self.edges.append({"data": data_dict})
        else:
            self.logger.debug("Move {} is already in moves".format(move))

    def jsonify_scope_viewer(self):
        """
        Use scope viewer to visualise pathways before the DBTL advances more.
        THe json file is a dict composed of one item called elements.
        The elements values is a dict composed of "nodes" and "edges"
        Nodes is a list of compounds, or reactions, with:
        """
        if self.pathway_as_dict is None:
            self.nodes_compounds.reverse()
            self.pathway_as_dict = {"elements": {"nodes": self.nodes_compounds + self.nodes_transformations,
                                            "edges": self.edges}}
        with open(self.file_to_save, "w") as json_handler:
            json.dump(self.pathway_as_dict, json_handler, indent = 2)

    def export_as_json_dict(self):
        """
        To export as a dict without needing to read and write the json.
        """
        if self.pathway_as_dict is None:
            self.nodes_compounds.reverse()
            self.pathway_as_dict = {"elements": {"nodes": self.nodes_compounds + self.nodes_transformations,
                                            "edges": self.edges}}
        return(self.pathway_as_dict)


def __cli():
    """Command line interface. Was actually used to make quick
    tests before implementing them in the testing file"""
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )
    logging.warning("CLI is not available for Pathway")


if __name__ == "__main__":
    __cli()
