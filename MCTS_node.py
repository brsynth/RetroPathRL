"""
Defines the Node of the MCTS:
- it is the basic object of a MCTS
- its state is a chemical compound state object
"""

# General utilities
import logging
import sys
import pickle

# RP3 specific classes
from compound import Compound, unpickle
from move import Move
from chemical_compounds_state import ChemicalCompoundState
from rewarding import Basic_Rollout_Reward
from UCT_policies import Biochemical_UCT_1, Nature_UCT, Classical_UCT_RAVE, Classical_UCT_with_bias, Classical_UCT
from Rollout_policies import Rollout_policy_first
from organisms import detectable_cmpds_H, ecoli_chassis_H, Test_organism_H
from organisms import detectable_cmpds_noH, ecoli_chassis_noH, Test_organism_noH

# General configuration
from config import *

class MCTS_node(object):
    """
    Defines what a node is:
    - has a state (set of chemical compounds) attribute at that stage of the retrosynthesis
    - state has the rules attribute
    - keeps the move and move score that generated it, as well as its parent
    - its level (depth). Currently used for visualisation but could be used for rewarding as well
    - the rewarding scheme considered when chassis is reached
    """

    logger = logging.getLogger(__name__)

    def __init__(self,
                 state, # state of Game
                 move = None,
                 history = None,
                 expansion_width = 5,
                 parent = None,
                 level = 0, # Generation information
                 maximum_depth = 10,  # Maximum allowed depth along the tree.
                 rewarding = Basic_Rollout_Reward,
                 progressive_bias_strategy = None,
                 main_layer_tree = True,
                 main_layer_chassis = True,
                 chemical_scorer = "RandomChemicalScorer",
                 biological_scorer = "RandomBiologicalScorer",
                 chemical_score = False,
                 virtual_visits = 0,
                 minimal_visit_counts = 1,
                 use_RAVE = False):
        self.move = move  # is a move obejct so score is encoded inside
        self.parent = parent  # Parent of the Node
        self.level = level
        self.maximum_depth = maximum_depth
        self.expansion_width = expansion_width
        self.chemical_score = chemical_score
        self.average_score = 0
        self.total_score = 0
        self.visits = virtual_visits
        self.virtual_visits = virtual_visits  # For child transmission
        self.minimal_visit_counts = minimal_visit_counts
        self.use_RAVE = use_RAVE
        self.organism = state.organism
        self.state = state.clone()
        self.state.set_chemical_scorer(chemical_scorer)
        self.state.set_biological_scorer(biological_scorer)
        self.chemical_scorer = chemical_scorer
        self.biological_scorer = biological_scorer

        self.root = (parent is None)  # Verification with other criterai below
        self.rewarding = rewarding
        # Chemical information:
        self.main_layer_tree = main_layer_tree
        self.main_layer_chassis = main_layer_chassis

        self.flag_for_extension = False

        # Check children from transposition table, not from here
        if use_transpositions:
            """
            Under development, should not be used.
            """
            self.set_hash_node()
            if self.hash in transposition_table.keys():
                reference_node = transposition_table[self.hash][0]
                self.children = reference_node.children
                self._helper_history(history)  # Update current history with information from father node
                self._helper_history(reference_node.history) # Update current history with information from reference_node
                reference_node._helper_history(self.history) # Update reference node history with merged information
                self.moves = reference_node.moves
                self.current_reward = reference_node.current_reward
                self.terminal = reference_node.terminal
                if reference_node.has_a_solved_child:
                    self.logger.warning("Merging nodes, from which one is solved, for node {}.".format(self.hash))
                self.has_a_solved_child = reference_node.has_a_solved_child
                self.total_score = reference_node.total_score
                self.average_score = 0
                transposition_table[self.hash].append(self)
            else:
                self.children = []
                self._helper_history(history)
                # Obtaining available moves
                self.moves = self.state.GetMoves(top_x = self.expansion_width, chemical_score = chemical_score)
                if retrosynthesis:
                    self.current_reward = self.state.GetResults_from_InChI_Keys(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
                elif biosensor:
                    self.current_reward = self.state.GetResultsForBiosensors(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
                else:
                    raise NotImplementedError
                no_available_moves = (self.moves == [])
                fully_chassis = (self.current_reward == self.rewarding.full_state_reward)
                self.terminal = (no_available_moves or fully_chassis or self.level == maximum_depth)
                self.has_a_solved_child = False
                transposition_table[self.hash] = [self]
        else:
            self.children = []
            # Defining a node's history
            # A node's history will be all compounds preivously generated that were not in the sink.
            self._helper_history(history)

            # Obtaining available moves
            self.moves = self.state.GetMoves(top_x = self.expansion_width, chemical_score = chemical_score)
            if retrosynthesis:
                self.current_reward = self.state.GetResults_from_InChI_Keys(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
            elif biosensor:
                self.current_reward = self.state.GetResultsForBiosensors(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
            else:
                raise NotImplementedError
            no_available_moves = (self.moves == [])
            fully_chassis = (self.current_reward == self.rewarding.full_state_reward)
            self.terminal = (no_available_moves or fully_chassis or self.level == maximum_depth)
            self.has_a_solved_child = False

        self.expanded_moves_ids = []
        self._progressive_bias(progressive_bias_strategy)
        self.progressive_bias_strategy = progressive_bias_strategy
        if use_toxicity:
            self.toxicity = min([cmp.toxicity for cmp in self.state.compound_list])

    def print_all_attributes(self):
        """Used for debugging: prints all attributes of interest"""
        print("Printing all attributes for node {}".format(self))
        print("Generating move {}".format(self.move))  # is a move obejct so score is ecnoded inside
        print("Parent {}".format(self.parent))  # Parent of the Node
        print("Current level is {} and max allowed depth {}".format(self.level, self.maximum_depth))
        print("expansion_width {}".format(self.expansion_width))
        print("children {}".format(self.children)) # Future node children
        print("average_score is {}, total is {} and visists are {}".format(self.average_score, self.total_score, self.visits))
        # print("Chassis organism is {}".format(self.organism))
        print("Current state is {}".format(self.state))
        print("Is root (y/n): {}".format(self.root)) # Verification with other criterai below
        print("Rewarding startegy {}".format(self.rewarding))
        # Chemical information:
        print("Main tree layer: {}, main chassis layer: {}".format(self.main_layer_tree, self.main_layer_chassis))
        print("Allowed moves here {}".format(self.moves))
        print("Is terminal (y/n) {}".format(self.terminal))
        print("Has a solved cild (y/n): {}".format(self.has_a_solved_child))


    def _helper_history(self, history):
        """
        Defines the node's history of previous states. Should be avoided in the search.
        """
        # self.history = []
        if history is None:
            history_list = []
        else:
            history_list = history.clone().compound_list
        for compound in self.state.compound_list:
            if self.organism is None:
                pass
            elif self.organism.compound_in_state(compound):
                pass
            else:
                if compound.in_list(history_list, main_layer = self.main_layer_chassis):
                    pass
                else:
                    history_list.append(compound)
        self.logger.debug("History list {} should at least contain a compound of the state {}".format(history_list, self.state))
        self.history = ChemicalCompoundState(history_list)

    def set_chemical_scorer(self, chemical_scorer):
        self.chemical_scorer = chemical_scorer

    def set_biological_scorer(self, biological_scorer):
        self.biological_scorer = biological_scorer

    def set_non_terminal(self):
        self.terminal = False

    def set_id(self, id):
        self.id = id  # This is for json viewer of tree

    def set_hash_node(self):
        """
        Not used at the moment.
        This function is used for transposition tables - for defining the key of the state.
        If we account for depth -ie: merged nodes need to have the same depth, we add it to the hash.
        """
        if use_transpositions_depth:
            self.hash = str(self.state) + '_' + str(self.level)
        else:
            self.hash = str(self.state)

    def __repr__(self):
        """ At the moment, returns its state. In the long run, can be changed to add counts"""
        return(str(self.state))

    def __str__(self):
        rep = ""
        if not self.move is None:
            rep = rep + "\t" * self.level + "Move name: {}/ Biological score : {}/ Chemical score: {}".format(self.move, round(self.move.biological_score, 2), round(self.move.chemical_score,2) ) + "\n"
        rep = rep + "\t" * self.level + "Count: {}/ Score : {}/ Children: {}/ Terminal : {}".format(self.visits, round(self.average_score, 2), len(self.children), self.terminal) + "\n"
        rep = rep + "\t" * self.level + str(self.state)
        return(rep)

    def __eq__(self, other):
        """
        Define equality between 2 Nodes
        2 Nodes are equal if:
        - same level (depth)
        - same creating move
        - same count visit
        - same total score
        - they have the same state
        """
        equal = (self.level == other.level)
        if not equal:
            return False
        if self.move is None:
            if other.move is None:
                pass
            else:
                self.logger.debug("One node has a generating move and not the other")
                return(False)
        else:
            if other.move is None:
                self.logger.debug("One node has a generating move and not the other")
            else:
                if self.main_layer_tree:
                    equal = (self.move.eq_main_layer(other.move))
                else:
                    equal = (self.move.eq_full_inchi_key(other.move))
        if not equal:
            return False
        equal = (self.visits == other.visits)
        if not equal:
            return False
        equal = (self.total_score == other.total_score)
        if not equal:
            return False
        equal = (self.total_score == other.total_score)
        if not equal:
            return False
        equal = (self.state == other.state)
        return(equal)

    def save(self, file_name = None, folder_address = "pickled_data"):
        if file_name is None:
            file_name = self.state.state_name + '-' + str(self.level) + '-' + str(self.visits)
        file_saving = open('{}/node_{}.pkl'.format(folder_address, file_name), 'wb')
        pickle.dump(self, file_saving)

    def _verif_initilisation(self):
        """
        Verify the Nodes are properly set up
        - root has no parent and one element
        - other nodes have parents and are not empty
        """
        is_proper_root = (self.parent is None and self.move is None and len(self.state) == 1 and self.level == 0)
        is_proper_node = (self.parent is not None and self.move is not None and len(self.state) >= 1 and self.level > 0)
        return(is_proper_root or is_proper_node)

    def _progressive_bias(self, progressive_bias_strategy):
        """
        Used if we want to implement progressive bias
        The idea would be to increase score on the states closest to the organism beforehand
        """
        if progressive_bias_strategy is None:
            self.progressive_bias = 0
        elif isinstance(progressive_bias_strategy, int):
            self.logger.info("Progressive bias with FPU of {}".format(progressive_bias_strategy))
            self.progressive_bias = progressive_bias_strategy
        elif progressive_bias_strategy == "state_reward":
            self.progressive_bias = self.current_reward
        elif progressive_bias_strategy == "max_reward":
            self.progressive_bias = self.rewarding.full_state_reward
        else:
            self.logger.error(type(progressive_bias_strategy))
            raise NotImplementedError

    def UCTSelectChild(self, UCT_policy):
        """
        Use UCT policy (defined in UCT policies) to select the best child.
        We are first ensuring all nodes have been visited at least minimal number of times.
        Then we use the UCT policy.
        """
        node_below_minimal = False
        if use_transpositions:
            self.children = transposition_table[self.hash][0].children
        for child in self.children:
            if child.visits <= self.minimal_visit_counts:
                node_below_minimal = True
                return(child)

        if not node_below_minimal:
            s = UCT_policy.calculate(self, top_n = 1)
            return(s)

    def SelectBestChild(self, policy = "visits"):
        """
        At the moment, the best child is the one with the highest number of visits.
        """
        if use_transpositions:
            self.children = transposition_table[self.hash][0].children
        if policy == "visits":
            s = sorted(self.children, key = lambda c: c.visits)
            s = s[-1]
            return(s)
        else:
            raise NotImplementedError

    def AddChild(self, move, state, remove = True):
        """
        Add a new child node for this move.
        If using transposition table, add it to the reference parent.
        """
        self.logger.info("Adding a child with move {} to node \n{}".format(move, self))
        state = state.clone()
        state.ApplyMove(move)
        # If untried moves list etc
        self.expanded_moves_ids.append(move)  # This is only useful if I want to extend tree at later iterations
        if remove:
            if use_transpositions:
                for node in transposition_table[self.hash]:
                    try:
                        node.moves.remove(move)
                    except ValueError:
                        self.logger.debug("For some reason, move {} already removed from node {}".format(move, node))
            else:
                self.moves.remove(move)  # Only if I go for the untried moves things
        if use_transpositions:
            n = MCTS_node(move = move,
                      parent = transposition_table[self.hash][0],
                      state = state.clone(),
                      level = self.level + 1,
                      rewarding = self.rewarding,
                      expansion_width = self.expansion_width,
                      maximum_depth = self.maximum_depth,
                      main_layer_tree = self.main_layer_tree,
                      main_layer_chassis = self.main_layer_chassis,
                      history = self.history,
                      chemical_scorer = self.chemical_scorer,
                      biological_scorer = self.biological_scorer,
                      chemical_score = self.chemical_score,
                      virtual_visits = self.virtual_visits,
                      progressive_bias_strategy = self.progressive_bias_strategy,
                      minimal_visit_counts = self.minimal_visit_counts,
                      use_rave = self.use_RAVE)
        else:
            n = MCTS_node(move = move,
                      parent = self,
                      state = state.clone(),
                      level = self.level + 1,
                      rewarding = self.rewarding,
                      expansion_width = self.expansion_width,
                      maximum_depth = self.maximum_depth,
                      main_layer_tree = self.main_layer_tree,
                      main_layer_chassis = self.main_layer_chassis,
                      history = self.history,
                      chemical_scorer = self.chemical_scorer,
                      biological_scorer = self.biological_scorer,
                      chemical_score = self.chemical_score,
                      virtual_visits = self.virtual_visits,
                      progressive_bias_strategy = self.progressive_bias_strategy,
                      minimal_visit_counts = self.minimal_visit_counts,
                      use_RAVE = self.use_RAVE)
        if use_transpositions:
            for node in transposition_table[self.hash]:
                node.children.append(n)
        else:
            self.children.append(n)
        if not self.terminal:
            self.terminal = (self.children == [])

    def expand(self):
        """
        Expand all children.
        Does not allow for moves that produce compounds that were present in the history
        """
        for move in self.moves:
            looping_move = False
            for product in move.product_list:
                if product.in_list(self.history.compound_list , main_layer = self.main_layer_chassis):
                    looping_move = True
                    self.logger.info("Looping move in node {}-{} for move {}".format(self.state, self.level, move))
                    break
            if not looping_move:
                self.AddChild(move, self.state, remove = False)
        if not self.terminal:
            self.terminal = (self.children == [])

    def expand_without_history(self):
        """
        Deprecated, now history is systematically verified.
        Expand all children.
        Will not change, but the GetMoves will to see how I handle moves
        Does not accoutn for the compound's history (ie: not regenrating previously seen compounds)
        """
        for move in self.moves:
            self.AddChild(move, self.state)

    def rollout(self, RolloutPolicy = Rollout_policy_first(), rollout_number = 5):
        """
        Performs a rollout on the current node
        - with a maximum depth of rollout_number
        - cloning the state before changing it dynamically
        """
        # Initilisation of the rollout:
        current_rollout = 0
        state = self.state.clone() # initial state before rollout
        available_moves = state.GetRolloutMoves()
        if retrosynthesis:
            in_chassis = (state.GetResults_from_InChI_Keys(rewarding = self.rewarding, main_layer = self.main_layer_chassis) == self.rewarding.full_state_reward)
        elif biosensor:
            in_chassis = (state.GetResultsForBiosensors(rewarding = self.rewarding, main_layer = self.main_layer_chassis) == self.rewarding.full_state_reward)
        else:
            raise NotImplementedError

        while (current_rollout < rollout_number and available_moves != [] and not in_chassis):
            move = RolloutPolicy.select_best_move(available_moves)
            state.ApplyMove(move)
            state = state.clone()
            if retrosynthesis:
                state_result = state.GetResults_from_InChI_Keys(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
            elif biosensor:
                state_result = state.GetResultsForBiosensors(rewarding = self.rewarding, main_layer = self.main_layer_chassis)
            else:
                raise NotImplementedError

            in_chassis = (state_result == self.rewarding.full_state_reward)
            current_rollout = current_rollout + 1
            available_moves = state.GetRolloutMoves()
        return(state)

    def update(self, result, solved, visit_number = 1, iteration = -1):
        """
        Update this node
        - visit_number additional visits
            - useful if parallelisation gets implemented
            - also is feedbakc adjustment gets implemented
        - total score
        - has_a_solved_child
        """
        if use_transpositions:
            # Visits are independent as else could cause huge bias.
            for node in transposition_table[self.hash]:
                node.visits = node.visits + visit_number  # This could change if there is too mcuh bias
                node.total_score = node.total_score + result * visit_number
                if not node.has_a_solved_child and solved: # if the state alreadyhas a solved child, well, it has one
                    node.has_a_solved_child = solved
                    node.logger.info("Udapted solved Node {} for the first time at iteration {}".format(node, iteration))
                try:
                    node.average_score = node.total_score/node.visits
                except ZeroDivisionError:
                    self.logger.warning("Division by zero when updating node {}".format(node))
                    node.average_score = node.total_score

        else:
            self.visits = self.visits + visit_number
            self.total_score = self.total_score + result * visit_number
            if not self.has_a_solved_child and solved: # if the state alreadyhas a solved child, well, it has one
                self.has_a_solved_child = solved
                self.logger.info("Udapted solved Node {} for the first time at iteration {}".format(self, iteration))
            self.average_score = self.total_score/self.visits
        if self.use_RAVE:
            if not self.move is None:
                self.move.update(result, visit_number)

    def flag_node_for_extension(self, extension_length = 10, maximum_depth = 15, chemical_scoring_configuration = None):
        """
        Used when extending from a previous search
        """
        self.flag_for_extension = True
        self.flag_extension_length = extension_length
        self.flag_maximum_depth = maximum_depth
        self.flag_chemical_scoring_configuration = chemical_scoring_configuration
        # Also reinitialise scores, because changing allowed rules and moves also changed status
        self.average_score = 0
        self.total_score = 0
        self.visits = self.virtual_visits

    def expand_later_iteration(self, extension_length = 10, maximum_depth = 15, chemical_scoring_configuration = None):
        """
        Deprecated.
        Was used when all nodes of the tree were extended
        Now nodes are flagged for extension and extended when visited, not at tree laoding.
        """

        self.maximum_depth = maximum_depth
        if self.level == self.maximum_depth:
            self.logger.warning("Did not extend node {} because it was a maximum depth of {}".format(self, self.level))
        self.expansion_width = self.expansion_width + extension_length
        for compound in self.state.compound_list:
            compound.set_chemical_scoring_configuration(chemical_scoring_configuration)
        all_moves = self.state.GetMoves(top_x = self.expansion_width, chemical_score = self.chemical_score, extension = True)
        for move in all_moves:
            if not move.in_list(self.expanded_moves_ids, self.main_layer_tree):
                self.moves.append(move)
        # self.logger.info("Extended node {} by {} moves".format(self.state, len(self.moves)))
        # Change terminal status
        no_available_moves = (self.moves == [])
        fully_chassis = (self.current_reward == self.rewarding.full_state_reward)
        if self.terminal:  # Only recalculate terminality for nodes that were considered terminal before
            self.terminal = (no_available_moves or fully_chassis or self.level == maximum_depth)
        self.flag_for_extension = False
        del self.flag_extension_length
        del self.flag_maximum_depth
        del self.flag_chemical_scoring_configuration
        self.logger.info("Extended node {} by {} moves. Terminality: {}".format(self.state, len(self.moves), self.terminal))


def __cli():
    """
    Command line interface.
    """
    print("There is no CLI available for MCTS objects")


if __name__ == "__main__":
    logging.basicConfig(
            stream=sys.stderr, level=logging.INFO,
            datefmt='%d/%m/%Y %H:%M:%S',
            format='%(asctime)s -- %(levelname)s -- %(message)s'
            )
    __cli()
