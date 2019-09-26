"""
Find supplements to complete a Tree.
Read argparser for details of arguments.
Principle is to identify compounds needed to complete chemical states.
"""

# General utilities
import os
import sys
import time
import signal
import datetime
import logging
import argparse
import pickle
import json

import random

from Tree import Tree

def unpickle(file_name, type = "tree", folder_address = "pickled_data"):
    with open('{}/{}_{}.pkl'.format(folder_address, type, file_name), 'rb') as input:
        return(pickle.load(input))


def run(tree, number_suggestions, rescued_states, folder_to_save, database = None):
    potential_supplements = {}
    # Extracting all potential supplements from the Tree.
    nodes_to_treat = [tree.root_node]
    while nodes_to_treat != []:
        node = nodes_to_treat[0]
        del nodes_to_treat[0]
        state = node.state
        supplement = state.GetSupplement_from_InChI_Keys()
        if not supplement is None:
            if supplement.InChIKey in potential_supplements.keys():
                potential_supplements[supplement.InChIKey]["rescued_states"] = potential_supplements[supplement.InChIKey]["rescued_states"] + 1
            else:
                information_to_keep = {"structure": supplement.csmiles,
                                       "name_from_MCTS": supplement.name,
                                       "synonyms_names": supplement.synonyms_names,
                                       "rescued_states":1}
                potential_supplements[supplement.InChIKey] = information_to_keep
        if node.terminal:
            pass
        else:
            for child in node.children:
                nodes_to_treat.append(child)
    logging.info("Potential supplements without filtering: {}".format(len(potential_supplements.keys())))
    # Sorting according to number of rescued states
    sorted_supplements = [suppl for suppl, value in sorted(potential_supplements.items(), key=lambda item: item[1]["rescued_states"], reverse=True) if value["rescued_states"] >= rescued_states]
    logging.info("Potential supplements after filtering with {} rescued states: {}".format(rescued_states, len(sorted_supplements)))

    # Filtering according to presence in a database of interest
    if database is None:
        supplements_of_interest = sorted_supplements
        logging.warning("Not checking availability within a Database of interest")
    else:
        supplements_of_interest = []
        for element in sorted_supplements:
            if element in database.keys():
                logging.info("Element {} (with {} pathways) is in database ({})".format(element, potential_supplements[element], database[element]))
                supplements_of_interest.append(element)
    # Filtering accoridng to maximal number of allwoed suggestions
    if len(supplements_of_interest) > number_suggestions:
        supplements_of_interest = supplements_of_interest[0:number_suggestions]
        logging.info("Keeping {} potential supplements".format(number_suggestions))
        assert len(supplements_of_interest) == number_suggestions
    else:
        logging.info("Keeping all supplements as there are only {} ({} allowed)".format(len(supplements_of_interest), number_suggestions))

    # Extracting pathways
    for supplement_to_extract in supplements_of_interest:
        # setting up search
        found_pathways = 0
        folder_to_save_pathways = "{}/{}".format(folder_to_save, supplement_to_extract.split("-")[0])
        if not os.path.exists(folder_to_save_pathways):
            os.mkdir(folder_to_save_pathways)
        # searching
        tree.set_folder_to_save(folder_to_save_pathways)
        nodes_to_treat = [tree.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            state = node.state
            supplement = state.GetSupplement_from_InChI_Keys()
            if not supplement is None:
                if supplement.InChIKey == supplement_to_extract:
                    found_pathways = found_pathways + 1
                    found_pathway = tree.extract_pathway_from_bottom(node, iteration=found_pathways)
            if node.terminal:
                pass
            else:
                for child in node.children:
                    nodes_to_treat.append(child)
        logging.info("Extract {} pathways for {}".format(found_pathways, supplement_to_extract))

def __cli():
    """
    Command line interface.
    """

    d = "Arguments for supplement finder. Find compounds that can complete a Tree and be suppelmented to media."
    parser = argparse.ArgumentParser(description=d)
    parser.add_argument("--tree_to_complete", help="Tree to find supplements to", default="end_search")
    parser.add_argument("--folder_tree_to_complete", help="Tree to find supplements to", default=None)

    parser.add_argument("--number_suggestions", default = 20,
                        help = "Maximum number of suggestions returned")
    parser.add_argument("--rescued_states", default = 1,
                        help = "Minimum number of times the compound must complete states")
    parser.add_argument("--folder_to_save", default="testing_supplement_finder")
    parser.add_argument("--terminal", help="Default logger is within the new folder_to_save, switch to terminal if specified",
                        action='store_true', default=False)
    parser.add_argument("--database_address", default=None,
                        help = "Address of a database to check availability. Json format required. Keys are inchikeys. Values are names, but could be cost or any metric of interest")

    args = parser.parse_args()
    folder_to_save = args.folder_to_save
    if not os.path.exists(folder_to_save):
        os.makedirs(folder_to_save, exist_ok=True)

    if args.terminal is True:
        logging.basicConfig(
                stream = sys.stderr,
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
    else:
        logging.basicConfig(
                stream = open("{}/{}.log".format(folder_to_save, "supplement_finder"), "w"),
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
    completed_tree = unpickle(file_name=args.tree_to_complete,
                           type='tree',
                        folder_address="{}/pickles".format(args.folder_tree_to_complete))
    if args.database_address is None:
        database = None
    else:
        with open(args.database_address, "r") as json_file:
            database = json.load(json_file)

    run(completed_tree, number_suggestions = args.number_suggestions,
        rescued_states =args.rescued_states, folder_to_save = args.folder_to_save,
        database = database)


if __name__ == "__main__":
    __cli()
