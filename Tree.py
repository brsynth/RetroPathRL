"""
Defines the Tree of the MCTS: main object for RetroPath 3
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
import string
import numpy as np
import json
import csv

import random

# RP3 specific packages
from move import Move  # only necessary for tracking class below
from rdkit import Chem  # only necessary for tracking below
from compound import Compound, unpickle, CompoundDefinitionException, ChemConversionError
from chemical_compounds_state import ChemicalCompoundState
from MCTS_node import MCTS_node
from representation import Test_representation, Test_to_file
from UCT_policies import Biochemical_UCT_1, Nature_UCT, Classical_UCT_RAVE, Classical_UCT_with_bias, Classical_UCT, \
    Biochemical_UCT_1_with_RAVE, Biochemical_UCT_with_progressive_bias, Chemical_UCT_1, Biological_UCT_1, \
    Biochemical_UCT_with_toxicity
from Rollout_policies import Rollout_policy_first, Rollout_policy_random_uniform_on_biochemical_multiplication_score, Rollout_policy_random_uniform
from rewarding import Basic_Rollout_Reward, RolloutRewards
from rule_sets_examples import applicable_rules_mixed_dict, applicable_rules_10_dict
from rule_sets_similarity import get_rules_and_score, full_rules_forward_H, full_rules_retro_H, full_rules_forward_no_H, \
    full_rules_retro_no_H
from pathway import Pathway
from pathway_scoring import RandomPathwayScorer, constant_pathway_scoring, null_pathway_scoring, \
    biological_pathway_scoring, biochemical_pathway_scoring
from tree_viewer import Tree_viewer
from organisms import detectable_cmpds_H, ecoli_chassis_H, Test_organism_H, iJO1366_chassis_H
from organisms import detectable_cmpds_noH, ecoli_chassis_noH, Test_organism_noH, import_organism_from_csv, iJO1366_chassis_noH
# General Configuration
from config import *
if use_toxicity:
    from compound_scoring import toxicity_scorer

sys.setrecursionlimit(50000)


class RunModeError(Exception):
    """Class for Tree Mode exception."""

    def __init__(self, retrosynthesis, biosensor):
        self.message = "Choose between retrosynthesis ({}) and biosensor ({})".format(retrosynthesis, biosensor)
        self.reason = "RunModeError"


class IncorrectTreeLoading(Exception):
    """Class for Conflicts between trees when loading a tree for search extension."""

    def __init__(self, message):
        self.message = message
        self.reason = "IncorrectTreeLoading"


class myTimeoutException(Exception):
    """Class for Tree Timeout exception."""

    def __init__(self, iteration=-1):
        self.iteration = iteration
        self.time = time.time()
        self.message = "Timeout at time {}, at iteration {}".format(self.time, self.iteration)


class FirstResultException(Exception):
    """Class for Tree First Result Stop exception."""

    def __init__(self, file_to_save, name, time_spent, iteration=-1):
        self.message = "First result after {}s for compound {}".format(time_spent, name)
        self.time = time_spent
        self.iteration = iteration
        with open(file_to_save, "w") as writer:
            writer.write("{};{}".format(name, round(time_spent, 2)))


class CompoundInSink(Exception):
    """Class for raising exception if compound already in sink."""

    def __init__(self, folder_to_save, name):
        self.message = "Compound {} already in organism".format(name)
        self.reason = "CompoundInSink"
        with open("{}/in_sink".format(folder_to_save), "w") as results_file:
            pass

class InvalidSink(Exception):
    """Class for raising exception if sink is invalid."""

    def __init__(self):
        self.reason = "InvalidSink"
        self.message = "Invalid (empty) sink"

# def timeout_tree_search(signum, frame):
#     """
#     Used for raising myTimeoutException when the alloted time has been spent.
#     Not in use at the moment - another strategy was defined that waits until the end of the iteration.
#     """
#     raise myTimeoutException

"""
The folllwoing commented functions were used for parallelisation.
They are not available at the publication of this tool.
"""

# def kill(pool):
#     """Send SIGTERMs to kill all processes belonging to pool.
#
#     Will not work on Windows OS.
#     Necessary for killing RDKit processes that don't listen when on C++
#     """
#     # stop repopulating new child
#     pool._state = mp.pool.TERMINATE
#     pool._worker_handler._state = mp.pool.TERMINATE
#     for p in pool._pool:
#         os.kill(p.pid, signal.SIGKILL)
#     # .is_alive() will reap dead process
#     while any(p.is_alive() for p in pool._pool):
#         pass
#     # Get-lucky workaround: force releasing lock
#     try:
#         pool._inqueue._rlock.release()
#     except ValueError as e:
#         logging.error(e)
#     pool.terminate()

# def worker_init_array(state_array, kwargs_local):
#     # Using code from https://stackoverflow.com/questions/39322677/python-how-to-use-value-and-array-in-multiprocessing-pool
#     global shared_array, kwargs
#     shared_array = state_array
#     kwargs = kwargs_local


# def worker_Rollout_with_array(i):
#     """
#     Used for storing intermediate results when performing multiple rollouts in parallel.
#     Dropped at the moment due to the necessity of isolating RDkit processes in their own workers.
#     - define multiple workers, each for a Rollout policy
#     - they would need to not require arguments, no I'd have to duplicate the biochemical ones
#     - the worker only takes the idnex of the array it is going to work on
#     - the array itself (shared memory) and the kwargs will be defined globallly
#     """
#
#     # kwargs = kwargs["kwargs"]
#     RolloutPolicy_name = kwargs["RolloutPolicy_name"]
#     if RolloutPolicy_name == "Rollout_policy_first":
#         RolloutPolicy = Rollout_policy_first()
#     elif RolloutPolicy_name == "Rollout_policy_chemical_best":
#         RolloutPolicy = Rollout_policy_chemical_best()
#     elif RolloutPolicy_name == "Rollout_policy_biological_best":
#         RolloutPolicy = Rollout_policy_biological_best()
#     elif RolloutPolicy_name == "Rollout_policy_biochemical_addition_best":
#         RolloutPolicy = Rollout_policy_biochemical_addition_best()
#     elif RolloutPolicy_name == "Rollout_policy_biochemical_multiplication_best":
#         RolloutPolicy = Rollout_policy_biochemical_multiplication_best()
#     elif RolloutPolicy_name == "Rollout_policy_random_uniform":
#         RolloutPolicy = Rollout_policy_random_uniform()
#     elif RolloutPolicy_name == "Rollout_policy_random_uniform_on_chem_score":
#         RolloutPolicy = Rollout_policy_random_uniform_on_chem_score()
#     elif RolloutPolicy_name == "Rollout_policy_random_uniform_on_bio_score":
#         RolloutPolicy = Rollout_policy_random_uniform_on_bio_score()
#     elif RolloutPolicy_name == "Rollout_policy_random_uniform_on_biochemical_addition_score":
#         RolloutPolicy = Rollout_policy_random_uniform_on_biochemical_addition_score()
#     elif RolloutPolicy_name == "Rollout_policy_random_uniform_on_biochemical_multiplication_score":
#         RolloutPolicy = Rollout_policy_random_uniform_on_biochemical_multiplication_score()
#     else:
#         raise NotImplementedError
#     state = kwargs["state"]
#     rewarding = kwargs["rewarding"]
#     rollout_number = kwargs["rollout_number"]
#     main_layer = kwargs["main_layer"]
#     current_rollout = 0
#     state = state.clone()  # initial state before rollout
#     available_moves = state.GetRolloutMoves(top_x=1)
#     if retrosynthesis:
#         state_result = state.GetResults_from_InChI_Keys(rewarding=rewarding, main_layer=main_layer)
#     elif biosensor:
#         state_result = state.GetResultsForBiosensors(rewarding=rewarding, main_layer=main_layer)
#     else:
#         raise NotImplementedError
#
#     in_chassis = (state_result == rewarding.full_state_reward)
#     shared_array[i] = state_result
#     while (current_rollout < rollout_number and available_moves != {} and not in_chassis):
#         move = RolloutPolicy.select_best_move(available_moves=available_moves)
#         # move = list(available_moves.keys())[0] #Warning: Rollout Policy needs to be pickled
#         state.ApplyMove(move)
#         state = state.clone()
#         if retrosynthesis:
#             state_result = state.GetResults_from_InChI_Keys(rewarding=rewarding, main_layer=main_layer)
#         elif biosensor:
#             state_result = state.GetResultsForBiosensors(rewarding=rewarding, main_layer=main_layer)
#         else:
#             raise NotImplementedError
#         in_chassis = (state_result == rewarding.full_state_reward)
#         current_rollout = current_rollout + 1
#         available_moves = state.GetRolloutMoves(top_x=1)
#         shared_array[i] = state_result
#     return state_result


class Tree(object):
    """
    Defines a search on the Tree according to
    - defined policy
    - organism
    Useful methods include:
    - running the search (key method of this object)
    - extracting the best pathway
    - extracting all pathways
    - visualising the tree
    - extending the tree with more tolerant settings or simply more iterations.
    - tree summary
    Numerous utilities has been developped for exporting results.
    """

    logger = logging.getLogger(__name__)

    def __init__(self,
                 root_state,  # starting chemical compound state (contains the target Compound only)
                 itermax=2,  # Iteration budget
                 time_budget=20, # Time budget
                 expansion_width=10,  # Number of allwoed children
                 available_rules=None,  # Rule dataset that can be used
                 rewarding=Basic_Rollout_Reward,  # How to reward the chemical state after rollout
                 max_depth=5,  # Maximum number of reactions to be used
                 max_rollout=3,
                 number_rollout=1,  # if parallel, number of different rollouts performed. NOT USED
                 UCT_policy="Classical_UCT",
                 UCT_parameters={"UCTK": 2, "bias_k": 0},
                 Rollout_policy="Rollout_policy_first",
                 Rollout_parameters={"rollout_timeout": 0.2},
                 print_tree_information=True,
                 parallel=False,
                 main_layer_tree=False,
                 main_layer_chassis=True,
                 organism=ecoli_chassis_H,
                 chemical_scorer="RandomChemicalScorer",
                 biological_scorer="RandomBiologicalScorer",
                 folder_to_save="tests/generated_jsons",
                 virtual_visits=0,
                 progressive_bias_strategy=None,
                 progressive_widening=False,
                 heavy_saving=False,
                 minimal_visit_counts=1,
                 use_RAVE=False,
                 pathway_scoring="null_pathway_scoring"):
        # Chemical equality between nodes in the stree: only main layer or full inchi
        self.main_layer_tree = main_layer_tree
        self.main_layer_chassis = main_layer_chassis
        self.organism = organism
        self.organism.set_main_layer(main_layer_chassis)
        self.root_state = root_state
        self.root_state.set_organism(self.organism)
        self.root_state.set_main_layer(self.main_layer_tree)  # Chassis only used in the chassis
        # Defaut rules for pytest - small subset
        if available_rules is None:
            self.available_rules = applicable_rules_10_dict
        else:
            self.available_rules = available_rules
        # Define rules and scores settings for the root state.
        self.root_state.set_available_rules(self.available_rules)
        self.root_state.set_chemical_scorer(chemical_scorer)
        self.root_state.set_biological_scorer(biological_scorer)

        self.scorer_info = "Chemical scoring is {}\n. Biological scoring is {}\n".format(chemical_scorer,
                                                                                         biological_scorer)
        # Define initial node of the Tree with all previous information
        self.root_node = MCTS_node(root_state,
                                   rewarding=rewarding,
                                   expansion_width=expansion_width,
                                   maximum_depth=max_depth,
                                   main_layer_tree=self.main_layer_tree,
                                   main_layer_chassis=self.main_layer_chassis,
                                   chemical_scorer=chemical_scorer,
                                   biological_scorer=biological_scorer,
                                   virtual_visits=virtual_visits,
                                   progressive_bias_strategy=progressive_bias_strategy,
                                   minimal_visit_counts=minimal_visit_counts,
                                   use_RAVE=use_RAVE)
        # Defining attributes from init
        self.expansion_width = expansion_width
        self.progressive_widening = progressive_widening
        self.itermax = itermax
        self.time_budget = time_budget
        self.rewarding = rewarding
        self.UCT_policy = UCT_policy
        self.UCT_parameters = UCT_parameters
        self.max_depth = max_depth
        self.max_rollout = max_rollout
        self.Rollout_policy = Rollout_policy
        self.Rollout_parameters = Rollout_parameters
        self.pathway_scoring = pathway_scoring
        # Saving settings
        self.folder_to_save = folder_to_save
        self.folder_to_save_pickles = "{}/{}".format(folder_to_save, "pickles")
        self.heavy_saving = heavy_saving

        self.print_tree_information = True

        self.parallel = parallel
        if self.parallel:
            """
            Should not be used at the moment.
            Used to conduct multiple rollouts in parallel.
            Impossible with current RDKit lack of timeout
            """
            self.worker = worker_Rollout_with_array
            self.number_rollout = number_rollout
        else:
            self.number_rollout = 1  # When not doing in parallel
        self.recap = None  # Used for storing results and information that will be outputed in simplified foramt.

    def set_folder_to_save(self, folder_to_save):
        self.folder_to_save = folder_to_save
        self.folder_to_save_pickles = "{}/{}".format(folder_to_save, "pickles")

    def set_heavy_saving(self, heavy_saving):
        self.heavy_saving = heavy_saving

    def set_rules(self, rules):
        self.available_rules = rules
        # Change also for all nodes.
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            node.state.set_available_rules(rules)
            for child in node.children:
                nodes_to_treat.append(child)

    def __len__(self):
        length = 0
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            length = length + 1
            for child in node.children:
                nodes_to_treat.append(child)
        return length

    def set_itermax(self, itermax):
        self.itermax = itermax

    def run_search(self, stop_at_first_result=False):
        """
        Main method for using RetroPath 3.
        Runs a tree search for itermax iterations:
        - starts from root
        - while node have children, select the best
        - expand the leaf if possible. Select the best child of the expanded node
        - rollout from this node
        - update all results
        - save pathways as they are found
        """
        UCT_policy = self._get_UCT_policy()
        pathway_scoring = self._get_pathway_scoring()
        if self.print_tree_information:
            self._print_tree_information_beginning()
        # Signal below is dropped in favour of checking at the beginning of an iteration if there is time left.
        # signal.signal(signal.SIGALRM, timeout_tree_search)
        # signal.alarm(self.time_budget)
        selection_time = 0
        expansion_time = 0
        rollout_time = 0
        update_time = 0
        raise_stop = False
        start_tree_search_time = time.time()
        try:
            for i in range(self.root_node.visits - self.root_node.virtual_visits,
                           self.itermax + self.root_node.visits - self.root_node.virtual_visits):
                # exploration and EXPANSION
                if time.time() - start_tree_search_time > self.time_budget:
                    raise myTimeoutException(iteration=i)
                self.logger.debug("Currently performing iteration {}".format(i))
                if i % int(round(self.itermax / 10)) == 0:
                    self.logger.info("Performing iteration {} out of {}".format(i, self.itermax))
                    if self.heavy_saving:
                        # Saving Tree for further analysis during teh search.
                        self.find_full_tree(iteration=i)
                        self.save(file_name="{}_iteration_{}".format(self.root_state.compound_list[0], i),
                                  folder_address=self.folder_to_save_pickles)
                        self.jsonify_full_tree(file_name="{}_iteration_{}".format(self.root_state.compound_list[0], i))
                # A new iteration starts from the root node.
                node = self.root_node
                start_selection_time = time.time()
                while node.children != [] and not node.terminal:
                    if node.flag_for_extension:
                        # When extending a previous tree.
                        node.expand_later_iteration(extension_length = node.flag_extension_length,
                                            maximum_depth = node.flag_maximum_depth,
                                            chemical_scoring_configuration = node.flag_chemical_scoring_configuration)

                    if not self.progressive_widening and node.moves == []:
                        # Go down the tree while all moves from a node has been expanded.
                        node = node.UCTSelectChild(UCT_policy=UCT_policy)
                        state = node.state
                    elif self.progressive_widening:
                        if node.moves == []: # No moves are available for PW.
                            node = node.UCTSelectChild(UCT_policy=UCT_policy)
                            state = node.state
                        else:  # there are still available moves, apply PW condition
                            if node.visits <= len(node.children) ^ 2:
                                node = node.UCTSelectChild(UCT_policy=UCT_policy)
                                state = node.state
                            else:
                                break
                    else:
                        break
                selection_time = selection_time + time.time() - start_selection_time
                start_expansion_time = time.time()

                if node.moves != [] and not node.terminal:  # This node can be expanded.
                    move = node.moves[0]
                    looping_move = False
                    for product in move.product_list:
                        if product.in_list(node.history.compound_list, main_layer=node.main_layer_chassis):
                            looping_move = True
                            self.logger.info("Looping move in node {}-{} for move {}".format(node.state, node.level, move))
                            node.moves.remove(move)
                            break
                    if not looping_move:
                        node.AddChild(move, node.state)
                        node = node.UCTSelectChild(UCT_policy=UCT_policy)
                        state = node.state

                expansion_time = expansion_time + time.time() - start_expansion_time
                # Perform rollout on the selected node
                if not node.terminal:
                    start_rollout_time = time.time()
                    allowed_rollout = min(self.max_rollout, self.max_depth - node.level)
                    if self.parallel:
                        """
                        Not in use at the moment.
                        Made for performing multiple rollouts on the same node.
                        Was temporarily suspended due to the necessity of encapsulating Rdkit preocesses.
                        """
                        shared_array_for_rollout = mp.Array('d', self.number_rollout)
                        kwargs_local = {"state": node.state,
                                        "rewarding": self.rewarding,
                                        "rollout_number": allowed_rollout,
                                        "RolloutPolicy_name": self.Rollout_policy,
                                        "main_layer": self.main_layer_chassis}
                        _pool = mp.Pool(processes=4, initializer=worker_init_array,
                                        initargs=(shared_array_for_rollout, kwargs_local))
                        try:
                            _pool.map(worker_Rollout_with_array, range(self.number_rollout))
                        except mp.TimeoutError as e:
                            self.logger.warning("Time out in parllel workers for node \n{}".format(node))
                            kill(_pool)
                        result, solved = np.mean(shared_array_for_rollout[:]), self.rewarding.full_state_reward in shared_array_for_rollout[:]
                    else:
                        # Perform rollout without parallelisation.
                        RolloutPolicy = self._get_Rollout_policy()
                        state = node.rollout(rollout_number=allowed_rollout, RolloutPolicy=RolloutPolicy)

                        if retrosynthesis:
                            result = state.GetResults_from_InChI_Keys(rewarding=self.rewarding,
                                                                      main_layer=self.main_layer_chassis)
                            solved = state.GetResults_from_InChI_Keys(rewarding=self.rewarding,
                                                                      main_layer=self.main_layer_chassis) == self.rewarding.full_state_reward
                        elif biosensor:
                            result = state.GetResultsForBiosensors(rewarding=self.rewarding,
                                                                   main_layer=self.main_layer_chassis)
                            solved = state.GetResultsForBiosensors(rewarding=self.rewarding,
                                                                   main_layer=self.main_layer_chassis) == self.rewarding.full_state_reward
                        else:
                            raise NotImplementedError
                    self.logger.info("Rollout results are {}, {}".format(result, solved))
                    rollout_time = rollout_time + time.time() - start_rollout_time
                else:
                    # Node is already terminal.
                    start_rollout_time = time.time()
                    state = node.state
                    if retrosynthesis:
                        result = state.GetResults_from_InChI_Keys(rewarding=self.rewarding,
                                                                  main_layer=self.main_layer_chassis)
                        solved = state.GetResults_from_InChI_Keys(rewarding=self.rewarding,
                                                                  main_layer=self.main_layer_chassis) == self.rewarding.full_state_reward
                    elif biosensor:
                        result = state.GetResultsForBiosensors(rewarding=self.rewarding,
                                                               main_layer=self.main_layer_chassis)
                        solved = state.GetResultsForBiosensors(rewarding=self.rewarding,
                                                               main_layer=self.main_layer_chassis) == self.rewarding.full_state_reward
                    else:
                        raise NotImplementedError
                    if solved and not node.has_a_solved_child:
                        # Extract pathway if encountered from the first time and fully solved.
                        self.logger.info("Pathway is found at iteration {}".format(i))
                        found_pathway = self.extract_pathway_from_bottom(node, iteration=i)
                        if len(found_pathway) == 1:
                            name = str(self.root_state.compound_list[0]) + "iteration_{}".format(i)
                            if self.heavy_saving:
                                found_pathway[0]["pathway"].save(name, folder_address=self.folder_to_save_pickles)
                            result = pathway_scoring.calculate(found_pathway[0]["pathway"])
                        else:
                            best_path = None
                            for path in found_pathway:
                                name = str(self.root_state.compound_list[0]) + "iteration_{}_{}".format(i,
                                                                                                    path["number"])
                                if self.heavy_saving:
                                    path["pathway"].save(name, folder_address=self.folder_to_save_pickles)
                                result = pathway_scoring.calculate(path["pathway"])
                                if best_path is None:
                                    best_path = path
                                    best_result = result
                                else:
                                    if best_result < result:
                                        best_path = path
                                        best_result = result
                        if stop_at_first_result:
                            raise_stop = True

                # Backpropagate:
                start_update_time = time.time()
                while node is not None:  # backpropagate from the expanded node and work back to the root node
                    node.update(result,
                                solved=solved,
                                visit_number=self.number_rollout,
                                iteration=i)
                    node = node.parent
                update_time = update_time + time.time() - start_update_time
                if raise_stop:
                    time_spent = time.time() - start_tree_search_time
                    raise FirstResultException(
                        file_to_save="{}/{}_execution_time.txt".format(self.folder_to_save, self.root_state.state_name),
                        name=self.root_state.state_name,
                        time_spent=time_spent,
                        iteration=i)
        except myTimeoutException as e:
            self.logger.warning(
                "Timeout happened on Tree search after {} iterations for a budget of {}s".format(i, self.time_budget))

            self.recap = {"TIME_EXECUTION": round(e.time - start_tree_search_time, 2),
                          "STOP_REASON": "timeout",
                          "NUMBER_ITERATION": e.iteration}
        except FirstResultException as e:
            self.logger.warning(e.message)
            self.recap = {"TIME_EXECUTION": round(e.time, 2),
                          "STOP_REASON": "first_result",
                          "NUMBER_ITERATION": e.iteration}
        if self.print_tree_information:
            self._print_tree_information_end(selection_time, expansion_time, rollout_time, update_time)
        logging.info("At the end of this tree search, the cache contains {}".format(len(home_made_cache.keys())))
        if self.recap is None:
            self.recap = {"TIME_EXECUTION": round(time.time() - start_tree_search_time, 2),
                          "STOP_REASON": "iteration",
                          "NUMBER_ITERATION": i}

        self.logger.warning("RECAP TIME_EXECUTION={}".format(self.recap["TIME_EXECUTION"]))
        self.logger.warning("RECAP STOP_REASON={}".format(self.recap["STOP_REASON"]))
        self.logger.warning("RECAP NUMBER_ITERATION={}".format(self.recap["NUMBER_ITERATION"]))

        if self.heavy_saving:
            self.find_full_tree(iteration=i)
            self.save(file_name="{}_iteration_{}".format(self.root_state.compound_list[0], i),
                      folder_address=self.folder_to_save_pickles)
            self.jsonify_full_tree(file_name="{}_iteration_{}".format(self.root_state.compound_list[0], i))

    def _get_UCT_policy(self):
        if self.UCT_policy == "Classical_UCT":
            policy = Classical_UCT(self.UCT_parameters)
        elif self.UCT_policy == "Classical_UCT_with_bias":
            policy = Classical_UCT_with_bias(self.UCT_parameters)
        elif self.UCT_policy == "Classical_UCT_RAVE":
            policy = Classical_UCT_RAVE(parameters=self.UCT_parameters)
        elif self.UCT_policy == "Nature_UCT":
            policy = Nature_UCT(self.UCT_parameters)
        elif self.UCT_policy == "Biochemical_UCT_1":
            policy = Biochemical_UCT_1(self.UCT_parameters)
        elif self.UCT_policy == "Biochemical_UCT_with_progressive_bias":
            policy = Biochemical_UCT_with_progressive_bias(self.UCT_parameters)
        elif self.UCT_policy == "Biochemical_UCT_1_with_RAVE":
            policy = Biochemical_UCT_1_with_RAVE(self.UCT_parameters)
        elif self.UCT_policy == "Chemical_UCT_1":
            policy = Chemical_UCT_1(self.UCT_parameters)
        elif self.UCT_policy == "Biological_UCT_1":
            policy = Biological_UCT_1(self.UCT_parameters)
        elif self.UCT_policy == "Biochemical_UCT_with_toxicity":
            policy = Biochemical_UCT_with_toxicity(self.UCT_parameters)
        else:
            raise NotImplementedError(self.UCT_policy)
        return policy

    def _get_Rollout_policy(self):
        if self.Rollout_policy == "Rollout_policy_first":
            policy = Rollout_policy_first()
        elif self.Rollout_policy == "Rollout_policy_chemical_best":
            policy = Rollout_policy_chemical_best()
        elif self.Rollout_policy == "Rollout_policy_biological_best":
            policy = Rollout_policy_biological_best()
        elif self.Rollout_policy == "Rollout_policy_biochemical_addition_best":
            policy = Rollout_policy_biochemical_addition_best()
        elif self.Rollout_policy == "Rollout_policy_biochemical_multiplication_best":
            policy = Rollout_policy_biochemical_multiplication_best()
        elif self.Rollout_policy == "Rollout_policy_random_uniform":
            policy = Rollout_policy_random_uniform()
        elif self.Rollout_policy == "Rollout_policy_random_uniform_on_chem_score":
            policy = Rollout_policy_random_uniform_on_chem_score()
        elif self.Rollout_policy == "Rollout_policy_random_uniform_on_bio_score":
            policy = Rollout_policy_random_uniform_on_bio_score()
        elif self.Rollout_policy == "Rollout_policy_random_uniform_on_biochemical_addition_score":
            policy = Rollout_policy_random_uniform_on_biochemical_addition_score()
        elif self.Rollout_policy == "Rollout_policy_random_uniform_on_biochemical_multiplication_score":
            policy = Rollout_policy_random_uniform_on_biochemical_multiplication_score()
        else:
            raise NotImplementedError(self.RolloutPolicy)
        return policy

    def _get_pathway_scoring(self):
        self.pathway_scoring_logs = False
        if self.pathway_scoring == "RandomPathwayScorer":
            pathway_scoring = RandomPathwayScorer
        elif self.pathway_scoring == "constant_pathway_scoring":
            pathway_scoring = constant_pathway_scoring
        elif self.pathway_scoring == "biological_pathway_scoring":
            pathway_scoring = biological_pathway_scoring
        elif self.pathway_scoring == "null_pathway_scoring":
            pathway_scoring = null_pathway_scoring
        else:
            raise NotImplementedError(self.pathway_scoring)
        return pathway_scoring

    def _print_tree_information_beginning(self):
        test_block = "---------------Beginning of search informations------------------ \n"
        text_block = test_block + "UCT policy is : {} \n".format(str(self.UCT_policy))
        text_block = text_block + "Rollout policy is : {} \n".format(str(self.Rollout_policy))
        text_block = text_block + "Rewarding policy is : {} \n".format(str(self.rewarding))
        text_block = text_block + 'Maximum number of iterations {} \n'.format(self.itermax)
        start_time = time.localtime()
        self.starttime = time.time()
        text_block = text_block + 'Tree search started running at {}h {}min {}s \n'.format(start_time.tm_hour,
                                                                                           start_time.tm_min,
                                                                                           start_time.tm_sec)
        text_block = text_block + 'Your tree search is currently running ... \n \n'
        self.logger.info(text_block)

    def _print_tree_information_end(self, selection_time=0, expansion_time=0, rollout_time=0, update_time=0):
        end_time = time.localtime()
        test_block = "---------------End of search informations------------------ \n"
        text_block = test_block + 'Tree search finished running at {}h {}min {}s \n'.format(end_time.tm_hour,
                                                                                            end_time.tm_min,
                                                                                            end_time.tm_sec)
        text_block = text_block + 'Your tree search took {:.1e}s \n \n'.format(time.time() - round(self.starttime))
        text_block = text_block + 'Time was spent in selection: {:.1e}s, expansion: {:.1e}s, rollout: {:.1e}s, update: {:.1e}s \n'.format(
            selection_time, expansion_time, rollout_time, update_time)
        if self.root_node.has_a_solved_child:
            text_block = text_block + "CONGRATULATIONS, THATS A WIN"
        else:
            text_block = text_block + "You just lost the game"
        self.logger.info(text_block)

    def __repr__(self):
        """
        Represents the tree horizontally for small cases
        """
        rep = ''
        node_list = [self.root_node]
        while node_list != []:
            current_node = node_list[-1]
            del node_list[-1]
            rep = rep + str(current_node.level) + str(current_node) + "\n"
            for child in current_node.children:
                node_list.append(child)
        return rep

    def __eq__(self, other):
        """
        Verifying whether 2 Trees are identical. Will be useful for comparison of UCT methods
        - Trees are identical if all their nodes are recursively identical
        """
        self_pile = [self.root_node]
        other_pile = [other.root_node]
        equal = (self.root_node == other.root_node)
        if not equal:
            return False
        while self_pile != [] and other_pile != [] and equal:
            current_node_self = self_pile[-1]
            del self_pile[-1]
            self_pile = self_pile + current_node_self.children
            found = False
            for node in other_pile:
                if node == current_node_self:
                    current_node_other = node
                    found = True
            if not found:
                # Quits the loop with the return statement.
                return False
            other_pile.remove(current_node_other)
            equal = (current_node_other == current_node_self)  # it's been selected for this
            other_pile = other_pile + current_node_other.children
        return equal

    def jsonify_full_tree(self, file_name=None, folder_to_save=None):
        """
        - prints the whole tree. Useful to see where the search got lost.
        - uses the json format for node and moves
        - create a Tree viewer object
        """

        if file_name is None:
            file_name = self.root_state.state_name
        if folder_to_save is None:
            folder_to_save = self.folder_to_save
        file_to_save = folder_to_save + "/{}_full_tree_for_MCTS.json".format(file_name)

        #  First, assign IDs to all moves and compounds
        id_move = 0
        id_node = 0
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            if node.move is not None:
                node.move.set_id(id_move)
                id_move = id_move + 1
            node.set_id(id_node)
            id_node = id_node + 1
            for child in node.children:
                nodes_to_treat.append(child)

        tree_viewer = Tree_viewer(file_to_save=file_to_save)
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            tree_viewer.add_node(node)
            for child in node.children:
                nodes_to_treat.append(child)
        tree_viewer.jsonify_tree_viewer()

    def save(self, file_name=None, folder_address="pickled_data"):
        if self.__class__.__module__ == "__main__":
            self.logger.warning("Tree should not be in main for pickling")
        if file_name is None:
            base_name = self.root_state.state_name
            i = 1
            self.logger.info("Tree is saved at {}/tree_{}_{}.pkl".format(folder_address, base_name, i))
            while os.path.exists("{}/tree_{}_{}.pkl".format(folder_address, base_name, i)):
                i += 1
            file_name = base_name + '_' + str(i)
        file_saving = open('{}/tree_{}.pkl'.format(folder_address, file_name), 'wb')
        pickle.dump(self, file_saving)

    def equality_visited_states(self, other):
        """
        Verifying whether 2 Trees explored the same states. Will be useful for comparison of UCT methods
        - Trees are identical if all their nodes are recursively identical
        - contrary to identify function, checks states but not visit counts.
        """
        self_pile = [self.root_node]
        other_pile = [other.root_node]
        equal = (self.root_node.state == other.root_node.state)
        if not equal:
            return False
        while self_pile != [] and other_pile != [] and equal:
            current_node_self = self_pile[-1]
            del self_pile[-1]
            self_pile = self_pile + current_node_self.children
            found = False
            for node in other_pile:
                if node.state == current_node_self.state:
                    current_node_other = node
                    found = True
            if not found:
                return False
            other_pile.remove(current_node_other)
            equal = (current_node_other.state == current_node_self.state)
            other_pile = other_pile + current_node_other.children
        return equal

    def extract_pathway_from_bottom(self, node, name=None, iteration=-1):
        """
        Extract pathways from found state (and not from the top down)
        """
        bottom_up_pathway = Pathway(first_iteration=iteration, target=None,
                                    compounds=[], moves=[],
                                    main_layer=self.main_layer_chassis,
                                    organism=self.organism,
                                    edges=[],
                                    nodes_compounds=[],
                                    nodes_transformations=[])
        target = self.root_state.compound_list[0]
        bottom_up_pathway.add_compound(target, is_source=1)
        if name is None:
            name = str(target)
        file_to_save = self.folder_to_save + "/{}_iteration_{}.json".format(name, iteration)
        bottom_up_pathway.set_file_to_save(file_to_save)
        if use_transpositions:
            complete_pathways = []
            pathways_to_treat = [{"pathway": bottom_up_pathway, "number": "a", "node_list": [node]}]
            while pathways_to_treat != []:
                current_pathway_node = pathways_to_treat[0]
                del pathways_to_treat[0]
                pathway = current_pathway_node["pathway"].clone()
                number = current_pathway_node["number"]
                node_list = current_pathway_node["node_list"]
                node = node_list[-1]
                if not node.parent is None:
                    node_parents = transposition_table[node.parent.hash]
                    if len(node_parents) == 1:
                        node_list.append(node.parent)
                        new_pathway = {"pathway": pathway,
                                       "number": number,
                                       "node_list": node_list}
                        pathways_to_treat.append(new_pathway)
                    else:
                        for i in range(len(node_parents)):
                            parent = node_parents[i]
                            new_node_list = copy.deepcopy(node_list)
                            new_pathway = pathway.clone()
                            new_node_list.append(parent)
                            new_pathway = {"pathway": new_pathway,
                                           "number": number + '_' + list(string.ascii_lowercase)[i],
                                           "node_list": new_node_list}
                            pathways_to_treat.append(new_pathway)
                else:
                    complete_pathways.append(current_pathway_node)
            for pathway_as_dict in complete_pathways:
                pathway = pathway_as_dict["pathway"]
                number = pathway_as_dict["number"]
                node_list = pathway_as_dict["node_list"]
                node_list.reverse()
                for top_down_node in node_list[1:]:
                    pathway.add_reaction(top_down_node.move, depth=top_down_node.level)
                file_to_save = self.folder_to_save + "/{}_iteration_{}_{}.json".format(name, iteration, number)
                pathway.set_file_to_save(file_to_save)
                pathway.jsonify_scope_viewer()
            return complete_pathways
        else:
            nodelist = [node]
            while node.parent is not None:
                # Going up the tree to get the parents
                nodelist.append(node.parent)
                node = node.parent
            nodelist.reverse()
            for top_down_node in nodelist[1:]:
                bottom_up_pathway.add_reaction(top_down_node.move, depth=top_down_node.level)
            bottom_up_pathway.set_file_to_save(file_to_save)
            bottom_up_pathway.jsonify_scope_viewer()
            return [{"pathway": bottom_up_pathway}]

    def find_single_best_pathway(self, policy="visits", folder_to_save="pathways", name=None):
        """
        Currently, best according to number of visits.
        """

        best_pathway = Pathway(first_iteration=-1, target=None,
                               compounds=[], moves=[],
                               main_layer=self.main_layer_chassis,
                               organism=self.organism,
                               edges=[],
                               nodes_compounds=[],
                               nodes_transformations=[])
        target = self.root_state.compound_list[0]
        best_pathway.add_compound(target, is_source=1)
        if policy == "visits":
            self.logger.info(
                "According to Wikipedia, choosing the pathway with the most visists \n This could be replaced by the best score")
        if name is None:
            name = str(target)
        file_to_save = folder_to_save + "/{}_best.json".format(name)
        best_pathway.set_file_to_save(file_to_save)
        node = self.root_node
        while not node.terminal:
            # Could be modified with children cashing
            node = node.SelectBestChild(policy=policy)
            best_pathway.add_reaction(node.move, depth=node.level)
        best_pathway.jsonify_scope_viewer()

    def find_multiple_best_pathways(self, folder_to_save="pathways",
                                    name=None, return_result=False,
                                    return_pathways = False, biochem_sorted = False):
        """
        Returns all solved pathways during the search.
        They are sorted according to biochemical score.
        """
        # Have a pile of patwhays to test and add. Only add patwhays that are fully solved. Export them all
        pathways_to_print = []
        pathway_iteration = 1
        if self.root_node.has_a_solved_child:
            initial_pathway = Pathway(first_iteration=-1, target=None,
                                      compounds=[], moves=[],
                                      main_layer=self.main_layer_chassis,
                                      organism=self.organism,
                                      edges=[],
                                      nodes_compounds=[],
                                      nodes_transformations=[])

            target = self.root_state.compound_list[0]
            initial_pathway.add_compound(target, is_source=1)
            if name is None:
                name = str(target)
            # Start by treating the initla pathway
            max_number = 0
            pathways_to_treat = [{"pathway": initial_pathway, "node": self.root_node, "number": 1}]
            while pathways_to_treat != []:
                # Obtain and delete patwhay
                current_pathway_node = pathways_to_treat[0]
                del pathways_to_treat[0]
                node = current_pathway_node["node"]
                pathway = current_pathway_node["pathway"]
                number = current_pathway_node["number"]
                if node.move is not None:
                    pathway.add_reaction(node.move, depth=node.level)
                if node.terminal:
                    # If the node is terminal, add the creating move in the pwathway and stop there
                    if pathway in pathways_to_print:
                        self.logger.debug("Pathway is already present: \n {}".format(pathway))
                    else:
                        pathways_to_print.append(pathway)
                else:
                    # Extracting the nodes that have solved children
                    nodes_to_treat = []
                    for child in node.children:
                        if child.has_a_solved_child:
                            nodes_to_treat.append(child)
                    children_with_solved_num = len(nodes_to_treat)
                    for i in range(len(nodes_to_treat)):
                        # for i in range(len(nodes_to_treat)):
                        new_pathway = pathway.clone()
                        new_number = max_number + i + number
                        pathways_to_treat.append(
                            {"pathway": new_pathway, "node": nodes_to_treat[i], "number": new_number})
                    max_number = new_number + 1
                current_pile_length = len(pathways_to_treat)
            self.logger.info("Found {} pathways for {}".format(len(pathways_to_print), str(target)))
        else:
            pathways_to_print = []
            self.logger.info("No pathway was found for this compound")
        if biochem_sorted:
            sorted_pathways_couples = [(pathway, round(biochemical_pathway_scoring.calculate(pathway), 2)) for pathway in pathways_to_print]
            sorted_pathways = sorted(sorted_pathways_couples, key=lambda item: item[1], reverse=True)
            for j in range(len(sorted_pathways)):
                pathway = sorted_pathways[j][0]
                file_to_save = folder_to_save + "{}_{}.json".format(name, j + 1)
                pathway.set_file_to_save(file_to_save)
                pathway.jsonify_scope_viewer()
            sorted_pathways = [patwhay for patwhay, score in sorted_pathways]
        else:
            sorted_pathways = pathways_to_print
            for j in range(len(sorted_pathways)):
                pathway = sorted_pathways[j]
                file_to_save = folder_to_save + "/{}_{}.json".format(name, j + 1)
                pathway.set_file_to_save(file_to_save)
                pathway.jsonify_scope_viewer()
        if return_result:
            return len(sorted_pathways)
        if return_pathways:
            return sorted_pathways

    def find_full_scope(self, folder_to_save="pathways", name=None):
        """
        Returns the scope of the compound: all pathways leading to the chassis.
        """
        pile_to_treat = []
        pathways_to_print = []
        pathway_iteration = 1

        full_scope = Pathway(first_iteration=-1, target=None,
                             compounds=[], moves=[],
                             main_layer=self.main_layer_chassis,
                             organism=self.organism,
                             edges=[],
                             nodes_compounds=[],
                             nodes_transformations=[])

        target = self.root_state.compound_list[0]
        full_scope.add_compound(target, is_source=1)
        if name is None:
            name = str(target)
        max_number = 0
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            #  print(node.terminal)
            if node.move is not None:
                full_scope.add_reaction(node.move, depth=node.level)
            if node.terminal:
                pass
            else:
                for child in node.children:
                    if child.has_a_solved_child:
                        nodes_to_treat.append(child)
        file_to_save = folder_to_save + "/{}_full_scope.json".format(name)
        full_scope.set_file_to_save(file_to_save)
        full_scope.jsonify_scope_viewer()

    def find_full_tree(self, folder_to_save=None, name=None, iteration=0):
        """
        Returns all compounds and moves visited during the search. Can be a big file.
        """
        pile_to_treat = []
        pathways_to_print = []
        if folder_to_save is None:
            folder_to_save = self.folder_to_save

        full_tree = Pathway(first_iteration=iteration, target=None,
                            compounds=[], moves=[],
                            main_layer=self.main_layer_chassis,
                            organism=self.organism,
                            edges=[],
                            nodes_compounds=[],
                            nodes_transformations=[])

        target = self.root_state.compound_list[0]
        full_tree.add_compound(target, is_source=1)
        if name is None:
            name = str(target)
        max_number = 0
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            if node.move is not None:
                full_tree.add_reaction(node.move, depth=node.level)
            if node.terminal:
                pass
            else:
                for child in node.children:
                    nodes_to_treat.append(child)
        file_to_save = folder_to_save + "/{}_full_tree_{}.json".format(name, iteration)
        full_tree.set_file_to_save(file_to_save)
        full_tree.jsonify_scope_viewer()

    def find_all_scopes(self, folder_to_save="pathways", name=None):
        for depth in range(1, self.max_depth + 1):
            self.find_specific_depth_scope(folder_to_save=folder_to_save, name=name, depth=depth)

    def find_specific_depth_scope(self, folder_to_save="pathways", name=None, depth=3):
        """
        Finds the scope up to a specific depth (less than n chemical reactions)
        """
        pile_to_treat = []
        pathways_to_print = []
        pathway_iteration = 1

        full_scope = Pathway(first_iteration=-1, target=None,
                             compounds=[], moves=[],
                             main_layer=self.main_layer_chassis,
                             organism=self.organism,
                             edges=[],
                             nodes_compounds=[],
                             nodes_transformations=[])

        target = self.root_state.compound_list[0]
        full_scope.add_compound(target, is_source=1)
        if name is None:
            name = str(target)
        max_number = 0
        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            if node.move is not None:
                full_scope.add_reaction(node.move, depth=node.level)
            if node.terminal or node.level >= depth:
                pass
            else:
                for child in node.children:
                    if child.has_a_solved_child:
                        nodes_to_treat.append(child)
        file_to_save = folder_to_save + "/{}_depth_{}_scope.json".format(name, depth)
        full_scope.set_file_to_save(file_to_save)
        full_scope.jsonify_scope_viewer()

    def flag_nodes_for_extension(self, extension_length=10, maximum_depth=None, chemical_scoring_configuration = None):
        """
        Flag the nodes that will have to be extended from the Loaded Tree.
        Function used only where starting search from an already explored tree
        """
        total_len = len(self)
        logging.info("Length of tree is {}".format(total_len))
        flagged_nodes = 0
        if maximum_depth is None:
            maximum_depth = self.max_depth
        else:
            self.logger.info("Changing maximum depth from {} to {}".format(self.max_depth, maximum_depth))

        nodes_to_treat = [self.root_node]
        while nodes_to_treat != []:
            node = nodes_to_treat[0]
            del nodes_to_treat[0]
            flagged_nodes = flagged_nodes + 1
            for child in node.children:
                nodes_to_treat.append(child)
            node.flag_node_for_extension(extension_length=extension_length, maximum_depth=maximum_depth,
                                         chemical_scoring_configuration = chemical_scoring_configuration)
        assert flagged_nodes == total_len
        logging.info("Flagged {} nodes".format(flagged_nodes))


def __cli():
    """
    Command line interface.
    All options for running RetroPath 3
    """

    def get_progressive_bias(progressive_bias_strategy):
        try:
            progressive_bias_strategy = int(progressive_bias_strategy)
        except ValueError:
            progressive_bias_strategy = progressive_bias_strategy
            if progressive_bias_strategy == "None":
                progressive_bias_strategy = None
        return progressive_bias_strategy

    def define_folder_to_save(folder):
        if folder is None:
            folder_to_save = os.path.join('debugging_results', args.c_name)
        else:
            folder_to_save = folder
        if not os.path.exists(folder_to_save):
            os.makedirs(folder_to_save, exist_ok=True)
        if not os.path.exists(os.path.join(folder_to_save, 'pickles')):
            os.mkdir(os.path.join(folder_to_save, 'pickles'))
        return folder_to_save

    def get_representation(representation):
        if representation:
            representation = Test_representation
        else:
            representation = Test_to_file
        return representation

    def get_organism(biosensor, organism_name = "none", complementary_sink=None, add_Hs=True):
        """
        Imports sinks.
        - detectable compounds for biosensors
        - all avaialble sinks.
        - complementary sink for media supplementation
        """
        if biosensor:
            logging.info("Using detectable compounds as sink as running in biosensor mode.")
            if add_Hs:
                organism = detectable_cmpds_H
            else:
                detectable_cmpds_noH
        else:
            if organism_name == "none":
                if complementary_sink is None:
                    logging.warning("Need to specify a sink")
                    raise InvalidSink
                else:
                    organism = import_organism_from_csv(complementary_sink, add_Hs=add_Hs)
                    logging.info("File {} is the sink".format(complementary_sink))
            elif organism_name == "test":
                if add_Hs:
                    organism = Test_organism_H
                else:
                    organism = Test_organism_noH
            elif organism_name == "ecoli":
                if add_Hs:
                    organism = ecoli_chassis_H
                else:
                    organism = ecoli_chassis_noH
            elif organism_name == "core_ecoli":
                if add_Hs:
                    organism = core_ecoli_H
                else:
                    organism = core_ecoli_noH
            elif organism_name == "ijo1366":
                if add_Hs:
                    organism = iJO1366_chassis_H
                else:
                    organism = iJO1366_chassis_noH
            elif organism_name == "bsubtilis":
                if add_Hs:
                    organism = bsubtilis_H
                else:
                    organism = bsubtilis_noH
            else:
                logging.warning("This organism is not implemented yet: {}".format(organism_name))
                raise NotImplementedError
        if not complementary_sink is None and organism_name != "none":
            cmpds_to_add = import_organism_from_csv(complementary_sink, add_Hs=add_Hs)
            organism.merge_states(cmpds_to_add)
            logging.info("Add compounds from {} to the sink".format(complementary_sink))
        return(organism)

    d = "All arguments to run a nice MCTS"
    parser = argparse.ArgumentParser(description=d)
    # Logs and saving information
    parser.add_argument("--verbose", help="Default logger is INFO, switch to DEBUG is specified",
                        dest='verbose', action='store_true', default=False)
    parser.add_argument("--log_file", help="Default logger is stderr, switch to log_file if specified",
                        default=None)
    parser.add_argument("--folder_to_save",
                        help="Folder to store results. Default: temp",
                        default="temp")
    parser.add_argument("--heavy_saving",
                        help='If True, will save the tree each max_iteration/10',
                        default=False,
                        type=lambda x: (str(x).lower() == 'true'))
    parser.add_argument("--stop_at_first_result",
                        help='If True, will stop the first time it encounters a fully solved pathway',
                        default=False,
                        type=lambda x: (str(x).lower() == 'true'))
    # Compound information
    parser.add_argument("--c_name", help="Compound name. Defaults to None (InchiKey)",
                        default=None)
    # One of the next 2 arguments has to be specified to give a target structure.
    parser.add_argument("--c_smiles", help="Compound smiles", default = None)
    parser.add_argument("--c_inchi", help="Compound inchi", default = None)
    # Timeouts on rdkit processes
    parser.add_argument("--fire_timeout", help = "Time allowed for firing one rule on one substrate",
                        type = float, default = 1)
    parser.add_argument("--standardisation_timeout", help = "Time allowed for standardising results from one application",
                        type = float, default = 5)
    # Complementary sink: if we need to supplement the media or use another sink than the provided ones.
    parser.add_argument("--organism_name", default = "ecoli",
                        choices = ["none", "test", "ecoli", "core_ecoli", "bsubtilis", "ijo1366"])
    parser.add_argument("--complementary_sink",
                        help="address of a csv file containing compounds to add",
                        default=None)
    # Visualisation
    parser.add_argument("--representation",
                        help="If activated, uses colors for representation. Otherwise, gedit compatible",
                        action='store_true', default=False)
    #  General MCTS parameters
    parser.add_argument("--itermax", help="Maximum number of tree iterations", default=1000, type=int)
    parser.add_argument("--parallel",
                        help="Using rollout parallelisation. Default is False. Should not be used at the moment.",
                        type=lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--expansion_width", help="Maximum number of children", default=5, type=int)
    parser.add_argument("--time_budget", help="Time budget", default=300, type=int)
    parser.add_argument("--max_depth", help="Maximum depth of search", default=7, type=int)
    parser.add_argument("--minimal_visit_counts", default=1, type=int,
                        help="Minimal number of times a node has to be rolled out before his brothers can be expanded")
    # UCT parameters
    parser.add_argument("--UCT_policy",
                        help="UCT policy for the tree search.",
                        choices=["Classical_UCT", "Biological_UCT_1",
                                 "Classical_UCT_with_bias", "Classical_UCT_RAVE",
                                 'Biochemical_UCT_with_progressive_bias', 'Biochemical_UCT_1_with_RAVE',
                                 "Biochemical_UCT_1", "Nature_UCT", "Chemical_UCT_1",
                                 "Biological_UCT_1", "Biochemical_UCT_with_toxicity"],
                        default="Classical_UCT")
    parser.add_argument("--UCTK",
                        help="UCTK for exploration/exploitation", type=float,
                        default=2)
    parser.add_argument("--bias_k",
                        help="bias_k for exploration/exploitation", type=float,
                        default=1)
    parser.add_argument("--k_rave",
                        help="k_rave for weighting of RAVE/MCTS. Number of visits before MCTS/RAVE = 50%%", type=float,
                        default=50)
    parser.add_argument("--use_RAVE",
                        help="Use rave or not",
                        type=lambda x: (str(x).lower() == 'true'), default=False)
    # Rewarding
    parser.add_argument("--penalty",
                        help="penalty for fully unsolved state",
                        type=int, default=-1)
    parser.add_argument("--full_state_reward",
                        help="full_state_reward for fully solved state",
                        type=int, default=2)
    parser.add_argument("--pathway_scoring",
                        help="pathway scoring function",
                        choices=["RandomPathwayScorer", "constant_pathway_scoring",
                                 "null_pathway_scoring", "biological_pathway_scoring"],
                        default="constant_pathway_scoring")
    # Rollout parameters
    parser.add_argument("--Rollout_policy",
                        help="Rollout_policy for the tree search.",
                        choices=["Rollout_policy_chemical_best",
                                 "Rollout_policy_random_uniform_on_biochemical_multiplication_score",
                                 "Rollout_policy_biological_best", "Rollout_policy_biochemical_addition_best",
                                 "Rollout_policy_biochemical_multiplication_best", "Rollout_policy_random_uniform",
                                 "Rollout_policy_random_uniform_on_chem_score",
                                 "Rollout_policy_random_uniform_on_bio_score",
                                 "Rollout_policy_random_uniform_on_biochemical_addition_score", "Rollout_policy_first"],
                        default="Rollout_policy_first")
    parser.add_argument("--max_rollout",
                        help="Max rollout number", type=int,
                        default=3)

    # Chemical and biological scoring
    parser.add_argument("--chemical_scoring",
                        help="Chemical scoring policy.",
                        choices=["RandomChemicalScorer",
                                 "SubstrateChemicalScorer",
                                 "SubandprodChemicalScorer",
                                 "ConstantChemicalScorer"],
                        default="SubandprodChemicalScorer")
    parser.add_argument("--biological_score_cut_off", default=0.1, type=float)
    parser.add_argument("--substrate_only_score_cut_off", default=0.3, type=float)
    parser.add_argument("--chemical_score_cut_off", default=0.3, type=float)
    # Bias parameters
    parser.add_argument("--virtual_visits",
                        help="Virtual visits", type=int,
                        default=0)
    parser.add_argument("--progressive_bias_strategy",
                        help="Progressive bias strategy",
                        default="max_reward")
    parser.add_argument("--progressive_widening",
                        help="progressive_widening",
                        type=lambda x: (str(x).lower() == 'true'), default=False)
    parser.add_argument("--diameter", nargs='+',
                        help="Diameters to consider", default=[16], type=int)
    parser.add_argument("--EC_filter", nargs='+',
                        help="EC numbers to consider for rules", default=None, type=str)
    parser.add_argument("--small", help="Use only a small subset", type=lambda x: (str(x).lower() == 'true'),
                        default=False)
    parser.add_argument("--seed",
                        help="Seed", type=int,
                        default=None)
    # Load from a previously run tree
    parser.add_argument("--tree_to_complete", help="Tree to restart the search from", default=None)
    parser.add_argument("--folder_tree_to_complete", help="Tree to restart the search from", default=None)

    args = parser.parse_args()
    # Config folder where to save data
    folder_to_save = define_folder_to_save(args.folder_to_save)

    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    if args.log_file is None:
        logging.basicConfig(stream=sys.stderr,
                            level=logging_level,
                            datefmt='%d/%m/%Y %H:%M:%S',
                            format='%(asctime)s -- %(levelname)s -- %(message)s')
    else:
        if not "log" in args.log_file:
            log_file = "log_" + args.log_file
        else:
            log_file = args.log_file
        log_writer = open("{}/{}".format(folder_to_save, log_file), "w")
        logging.basicConfig(stream=log_writer,
                            level=logging_level,
                            datefmt='%d/%m/%Y %H:%M:%S',
                            format='%(asctime)s -- %(levelname)s -- %(message)s')
    if not args.seed is None:
        random.seed(args.seed)
        logging.warning("Setting the seed at {}".format(args.seed))

    # Information from argparse
    representation = get_representation(args.representation)
    rules, biological_scoring = get_rules_and_score(full_rules_forward_H=full_rules_forward_H,
                                                    full_rules_retro_H=full_rules_retro_H,
                                                    full_rules_forward_no_H=full_rules_forward_no_H,
                                                    full_rules_retro_no_H=full_rules_retro_no_H,
                                                    add_Hs=add_Hs,
                                                    retro=retrosynthesis,
                                                    diameters=args.diameter,
                                                    small=args.small,
                                                    c_name=args.c_name,
                                                    filtering_EC=args.EC_filter)
    if not args.EC_filter is None:
        logging.info("Filtering rules based on EC number. Currently {} rules".format(len(rules.keys())))

    progressive_bias_strategy = get_progressive_bias(args.progressive_bias_strategy)
    # Rewarding:
    rollout_rewards = RolloutRewards(penalty=args.penalty, full_state_reward=args.full_state_reward)
    # Stop mode
    if args.stop_at_first_result:
        stop_mode_config = "Stopping at first result"
    else:
        stop_mode_config = "Full search - no stopping at first result"
    # Minimal visits
    minimal_visit_counts_config = "Setting the minimal visit count for a node at {}".format(args.minimal_visit_counts)
    # RAVE_config
    RAVE_config = "Using RAVE: {}".format(args.use_RAVE)
    # Scoring configuration
    chemical_scoring = args.chemical_scoring
    chemical_scoring_configuration = {
        "biological_score_cut_off": args.biological_score_cut_off,
        "substrate_only_score_cut_off": args.substrate_only_score_cut_off,
        "chemical_score_cut_off": args.chemical_score_cut_off}

    biological_score_config = "Using biological cut off at {}".format(args.biological_score_cut_off)
    substrate_only_score_config = "Using substrate only cut off at {}".format(args.substrate_only_score_cut_off)
    chemical_score_config = "Using chemical score cut off at {}".format(args.chemical_score_cut_off)

    # Setting chemistry info
    logging.info("Stating global parameters from configuration file")
    logging.info(tree_mode_config)
    logging.info(stop_mode_config)
    logging.info(DB_config)
    logging.info(cache_config)
    logging.info("-------------------Scoring configurations -------------------------")
    logging.info(biological_score_config)
    logging.info(substrate_only_score_config)
    logging.info(chemical_score_config)
    logging.info("-------------------MCTS configurations -------------------------")
    logging.info(minimal_visit_counts_config)
    logging.info(RAVE_config)
    logging.info("-------------------Chemical configurations -------------------------")
    logging.info(hydrogen_config)
    if use_toxicity:
        logging.info("-------------------Toxicity configurations -------------------------")
        logging.info(toxicity_scorer.log_loading)
        logging.info(toxicity_scorer.log_score)
    try:
        organism = get_organism(biosensor, organism_name = args.organism_name, complementary_sink = args.complementary_sink, add_Hs = add_Hs)
        if retrosynthesis and biosensor:
            raise RunModeError(retrosynthesis, biosensor)
        try:
            if args.tree_to_complete is None:
                root_compound = Compound(csmiles = args.c_smiles,
                                        InChI = args.c_inchi,
                                        name = args.c_name,
                                        max_moves = args.expansion_width,
                                        stereo = False,
                                        heavy_standardisation = True,
                                        fire_timeout = args.fire_timeout,
                                        chemical_scoring_configuration = chemical_scoring_configuration,
                                        standardisation_timeout = args.standardisation_timeout)
                state = ChemicalCompoundState([root_compound], representation = representation)  # state is not sanitised
                if biosensor:
                    present_in_state_detectable = organism.compound_in_state(root_compound)
                    if present_in_state_detectable:
                        logging.warning("Removed compound from the detectable set to force enzymatic detection")
                        organism.remove_cmpd_from_state(root_compound)
                else:
                    present_in_state_sink = organism.compound_in_state(root_compound)
                    if present_in_state_sink:
                        raise CompoundInSink(folder_to_save, root_compound)

                search_tree = Tree(root_state=state,
                                   itermax=args.itermax,
                                   parallel=args.parallel,
                                   available_rules=rules,
                                   rewarding=rollout_rewards,
                                   expansion_width=args.expansion_width,
                                   time_budget=args.time_budget,
                                   max_depth=args.max_depth,
                                   UCT_policy=args.UCT_policy,
                                   UCT_parameters={"UCTK": args.UCTK, "bias_k": args.bias_k, 'k_rave': args.k_rave},
                                   Rollout_policy=args.Rollout_policy,
                                   max_rollout=args.max_rollout,
                                   organism=organism,
                                   main_layer_tree=True,
                                   main_layer_chassis=True,
                                   biological_scorer=biological_scoring,
                                   chemical_scorer=chemical_scoring,
                                   folder_to_save=folder_to_save,
                                   virtual_visits=args.virtual_visits,
                                   progressive_bias_strategy=progressive_bias_strategy,
                                   progressive_widening=args.progressive_widening,
                                   heavy_saving=args.heavy_saving,
                                   minimal_visit_counts=args.minimal_visit_counts,
                                   use_RAVE=args.use_RAVE,
                                   pathway_scoring=args.pathway_scoring)
            else:
                search_tree = unpickle(file_name=args.tree_to_complete,
                                       type='tree',
                                       folder_address="{}/pickles".format(args.folder_tree_to_complete))
                # Check compound compatibility
                current_root_state = search_tree.root_state
                try:
                    root_compound = Compound(csmiles=args.c_smiles,
                                             InChI=args.c_inchi,
                                             name=args.c_name,
                                             max_moves=args.expansion_width,
                                             stereo=False,
                                             heavy_standardisation=True)
                    state = ChemicalCompoundState([root_compound],
                                                  representation=representation)  # state is not sanitised
                    if state != current_root_state:
                        raise IncorrectTreeLoading(
                            "New root {} is different from old root {} when loading tree".format(root_compound,
                                                                                                 current_root_state.compound_list[
                                                                                                     0]))
                except CompoundDefinitionException:
                    logging.warning("Use compound information from previous Tree")
                    root_compound = current_root_state
                if biosensor:
                    present_in_state_detectable = organism.compound_in_state(root_compound)
                    if present_in_state_detectable:
                        logging.warning("Removed compound from the detectable set to force enzymatic detection")
                        organism.remove_cmpd_from_state(root_compound)

                search_tree.set_heavy_saving(args.heavy_saving)
                search_tree.set_folder_to_save(folder_to_save)
                search_tree.find_full_scope(folder_to_save=folder_to_save, name="after_unpickling")
                search_tree.jsonify_full_tree(file_name="after_unpickling")

                search_tree.set_rules(rules)  # Resetting the rules to new standards and not former tree
                start_time = time.time()
                logging.info("Starting flagging for extension at {}".format(start_time))
                search_tree.flag_nodes_for_extension(extension_length=args.expansion_width,
                                             maximum_depth=args.max_depth,
                                             chemical_scoring_configuration = chemical_scoring_configuration)
                logging.info("Finished flagging for extension at {}".format(time.time() - start_time))
                search_tree.find_full_scope(folder_to_save=folder_to_save, name="after_extending_nodes")
                search_tree.jsonify_full_tree(file_name="after_extending_nodes")
            # Running the search for both trees (new or loaded)
            search_tree.run_search(args.stop_at_first_result)
        except KeyboardInterrupt as e:
            logging.warning("Keyboard interruption")

        search_tree.jsonify_full_tree()
        search_tree.find_full_scope(folder_to_save=folder_to_save)
        if args.heavy_saving:
            search_tree.find_all_scopes(folder_to_save=folder_to_save)
        try:
            search_tree.find_single_best_pathway(folder_to_save=folder_to_save)
        except IndexError:
            logging.info("No best pathway was found")
        nbr = search_tree.find_multiple_best_pathways(folder_to_save=folder_to_save, return_result=True)
        logging.info("{} pathways were found".format(nbr))
        logging.info("At the end of this tree search, the cache contains {}".format(len(home_made_cache.keys())))
        search_tree.save(folder_address=folder_to_save + '/pickles', file_name="end_search")
    except (RunModeError, IncorrectTreeLoading, CompoundInSink, InvalidSink) as e:
        logging.error(e.message)
        nbr = 0
        loading_recap = {"TIME_EXECUTION": 0, "STOP_REASON": e.reason, "NUMBER_ITERATION": 0}
        logging.warning("RECAP TIME_EXECUTION={}".format(loading_recap["TIME_EXECUTION"]))
        logging.warning("RECAP STOP_REASON={}".format(loading_recap["STOP_REASON"]))
        logging.warning("RECAP NUMBER_ITERATION={}".format(loading_recap["NUMBER_ITERATION"]))
    except ChemConversionError as e:
        logging.error(e)
        nbr = 0
        loading_recap = {"TIME_EXECUTION": 0, "STOP_REASON": "incorrect_input", "NUMBER_ITERATION": 0}
        logging.error("Verify the input SMILES or InChI is valid")
        logging.warning("RECAP TIME_EXECUTION={}".format(loading_recap["TIME_EXECUTION"]))
        logging.warning("RECAP STOP_REASON={}".format(loading_recap["STOP_REASON"]))
        logging.warning("RECAP NUMBER_ITERATION={}".format(loading_recap["NUMBER_ITERATION"]))

    # Exporting information on obtained results
    with open("{}/results.csv".format(folder_to_save), "w") as csv_file:
        csv_writer = csv.DictWriter(csv_file, fieldnames=["parameter", "value"])
        csv_writer.writeheader()
        csv_writer.writerow({"parameter": "stop_at_first_result", "value": args.stop_at_first_result})
        csv_writer.writerow({"parameter": "c_name", "value": args.c_name})
        csv_writer.writerow({"parameter": "c_smiles", "value": args.c_smiles})
        csv_writer.writerow({"parameter": "c_inchi", "value": args.c_inchi})
        csv_writer.writerow({"parameter": "fire_timeout", "value": args.fire_timeout})
        csv_writer.writerow({"parameter": "organism_name", "value": args.organism_name})
        csv_writer.writerow({"parameter": "complementary_sink", "value": args.complementary_sink})
        csv_writer.writerow({"parameter": "itermax", "value": args.itermax})
        csv_writer.writerow({"parameter": "expansion_width", "value": args.expansion_width})
        csv_writer.writerow({"parameter": "time_budget", "value": args.time_budget})
        csv_writer.writerow({"parameter": "max_depth", "value": args.max_depth})
        csv_writer.writerow({"parameter": "minimal_visit_counts", "value": args.minimal_visit_counts})
        csv_writer.writerow({"parameter": "UCT_policy", "value": args.UCT_policy})
        csv_writer.writerow({"parameter": "UCTK", "value": args.UCTK})
        csv_writer.writerow({"parameter": "bias_k", "value": args.bias_k})
        csv_writer.writerow({"parameter": "k_rave", "value": args.k_rave})
        csv_writer.writerow({"parameter": "use_RAVE", "value": args.use_RAVE})
        csv_writer.writerow({"parameter": "penalty", "value": args.penalty})
        csv_writer.writerow({"parameter": "full_state_reward", "value": args.full_state_reward})
        csv_writer.writerow({"parameter": "Rollout_policy", "value": args.Rollout_policy})
        csv_writer.writerow({"parameter": "max_rollout", "value": args.max_rollout})
        csv_writer.writerow({"parameter": "chemical_scoring", "value": args.chemical_scoring})
        csv_writer.writerow({"parameter": "biological_score_cut_off", "value": args.biological_score_cut_off})
        csv_writer.writerow({"parameter": "substrate_only_score_cut_off", "value": args.substrate_only_score_cut_off})
        csv_writer.writerow({"parameter": "chemical_score_cut_off", "value": args.chemical_score_cut_off})
        csv_writer.writerow({"parameter": "virtual_visits", "value": args.virtual_visits})
        csv_writer.writerow({"parameter": "progressive_bias_strategy", "value": args.progressive_bias_strategy})
        csv_writer.writerow({"parameter": "progressive_widening", "value": args.progressive_widening})
        csv_writer.writerow({"parameter": "diameter", "value": args.diameter})
        csv_writer.writerow({"parameter": "EC_filter", "value": args.EC_filter})
        csv_writer.writerow({"parameter": "tree_to_complete", "value": args.tree_to_complete})
        csv_writer.writerow({"parameter": "found_pathways", "value": nbr})
        try:
            recap = search_tree.recap
        except UnboundLocalError:
            recap = loading_recap
        csv_writer.writerow({"parameter": "TIME_EXECUTION", "value": recap["TIME_EXECUTION"]})
        csv_writer.writerow({"parameter": "STOP_REASON", "value": recap["STOP_REASON"]})
        csv_writer.writerow({"parameter": "NUMBER_ITERATION", "value": recap["NUMBER_ITERATION"]})

    with open("{}/{}_results".format(folder_to_save, nbr), "w") as results_file:
        # Only printing out results numbers.
        pass


if __name__ == "__main__":
    from Tree import *
    __cli()
