"""
Defines the state of the MCTS:
It is the basic object inside a node. It contains all compounds necessary for retrosynthesis.
- set of compounds, and moves that apply to them.
- moves are available transformations.
- all transformations on all compounds are considered and only the best are used for expansion.
"""

# General utilities
import logging
import sys
import pickle
import bisect  # Faster insertion in ordered lists

# Object classes for RP3
from compound import Compound, unpickle  # Base chemical compound object
from representation import Test_representation, Test_to_file  # How to represent the results (ie: colors or not)
from rewarding import Basic_Rollout_Reward  # How rollout runs/ chassis compounds are rewarded
from rule_sets_examples import applicable_rules_10_dict  # For an example

# General configuration
from config import *

class ChemicalStateException(Exception):
    """Base for this object's exceptions.."""
    pass

class ChemicalSanitisationException(ChemicalStateException):
    """Raised when the state is not properly sanitised."""

    def __init__(self):
        self._msg = "CHEMICALSTATE-SANITISATION-ERROR"

    def __str__(self):
        return self._msg


class NotImplementedError(Exception):
    """Raised when a function is not yet implemented"""

    def __init__(self, msg = "Not implemented yet"):
        self._msg = msg

    def __str__(self):
        return self._msg

class ChemicalStateNamingException(ChemicalStateException):
    """Raised when there is an error in naming the state and its compound"""

    def __init__(self, msg = "Error while naming"):
        self._msg = msg

    def __str__(self):
        return self._msg


class ChemicalCompoundState(object):
    """
    Defines what a state is and how to perform moves.
    A state contains:
    - compounds (initialised by compound list, possibly named)
    - an organism for States that are not themselves organisms. This allows reward function calculations
    - representation: how to represent elements of the state
    - an optional state name (for future transformation hashing)
    - a method to obtain potential supplements
    """

    logger = logging.getLogger(__name__)

    def __init__(self,
                 compound_list,
                 name_list = None,
                 organism = None,
                 representation = Test_representation,
                 state_name = None,
                 available_rules = None,
                 main_layer = True,
                 automatic_sanitation = True,
                 chemical_scorer = "RandomChemicalScorer",
                 biological_scorer = "RandomBiologicalScorer"):

        self.main_layer = main_layer
        self._helper_init(compound_list, name_list) # Treats various ways to Initialize the state
        self.sanitised = self._check_sanitized_state()  # Verify it has been properly set up
        self.length = len(self.compound_list)
        self.organism = organism
        self.representation = representation
        self.chemical_scorer = chemical_scorer
        self.biological_scorer = biological_scorer
        if state_name is None:
            self._naming()
        else:
            # Usually not used
            self.state_name = state_name
        if available_rules is None:
            self.available_rules = applicable_rules_10_dict
        else:
            self.available_rules = available_rules
        if (not self.sanitised) and automatic_sanitation:
            self._sanitise()

    def _helper_init(self, compound_list, name_list):
        """
        Allows for set-up with:
        - Compound
        - Smile
        - List of compounds
        - list of smiles
        """
        if type(compound_list) is Compound:
            self.compound_list = [compound_list]
        elif type(compound_list) is list:
            if compound_list == []:
                self.compound_list = []
                self.logger.warning("Generating empty compound state")
            elif type(compound_list[0]) is Compound:
                self.compound_list = compound_list
            elif type(compound_list[0]) is str:
                self.compound_list = [Compound(x) for x in compound_list]
            else:
                print("Type of list elements is {}".format(type(compound_list[0])))
        elif type(compound_list) is str:
            self.compound_list = [Compound(compound_list)]
        else:
            print(type(compound_list))

        if name_list is not None:
            if type(name_list) is str:
                if len(self.compound_list) != 1:
                    raise ChemicalStateNamingException("Name list is a string and compound list is of len {}".format(len(self.compound_list)))
                else:
                    self.compound_list[0].naming(name_list)
            elif type(name_list) is list:
                if len(self.compound_list) != len(name_list):
                    raise ChemicalStateNamingException("Name and compound lists are of different length")
                else:
                    for i in range(len(compound_list)):
                        self.compound_list[i].naming(name_list[i])
        # sort the list
        self.compound_list = sorted(self.compound_list)

    def __len__(self):
        """Define state len as the number of compounds"""
        return len(self.compound_list)

    def __eq__(self, other):
        """
        Defines when 2 states are equal:
        - both sanitised
        - same length
        - all compounds of 1 are in the other
        - definitions of compounds IN the sink changes: main layer or full inchikey
        """
        # If states are not sanitised, no point in checking further.
        equal = self.sanitised and other.sanitised
        self.logger.debug("self.sanitised {}".format(self.sanitised))
        self.logger.debug("other.sanitised {}".format(other.sanitised))
        if not equal:
            return(equal)
        # If length are different, states are different
        # equal = equal and (len(self.compound_list) == len(other.compound_list))
        equal = len(self.compound_list) == len(other.compound_list)
        if not equal:
            self.logger.info("State lengths are not equal. Self's main layer is {} and other is {}".format(self.main_layer, other.main_layer))
            return(equal)
        # If all compounds of set one are in the other and length are equal, then they are equal
        for compound in self.compound_list:
            equal = equal and (compound.in_list(list_of_compounds = other.compound_list, main_layer = self.main_layer))
        return(equal)

    def set_available_rules(self, available_rules):
        self.available_rules = available_rules

    def set_chemical_scorer(self, chemical_scorer):
        self.chemical_scorer = chemical_scorer

    def set_biological_scorer(self, biological_scorer):
        self.biological_scorer = biological_scorer

    def set_main_layer(self, main_layer):
        """
        If True, use only the main layer from InChIKey, instead of the full InChIKey
        """
        self.main_layer = main_layer
        if self.main_layer:
            # If using only the main layer of the inchi, can remove duplicates.
            self.logger.info("Main layer was set to True. Sanitising state")
            self._sanitise()

    def set_organism(self, organism):
         self.organism = organism

    def save(self, file_name = None, folder_address = "pickled_data"):
        if file_name is None:
            file_name = self.state_name
        file_saving = open('{}/state_{}.pkl'.format(folder_address, file_name), 'wb')
        pickle.dump(self, file_saving)

    def not_in_self(self, other):
        """
        Returns all compounds that are in other state and not in self.
        """
        not_in_self = []
        for compound in other.compound_list:
            if not compound.in_list(self.compound_list, main_layer = self.main_layer):
                not_in_self.append(compound)
        return(not_in_self)

    def not_in_other(self, other):
        """
        Returns all compounds that are in self state and not in the other set.
        """
        not_in_other = []
        for compound in self.compound_list:
            if not compound.in_list(other.compound_list, main_layer = self.main_layer):
                not_in_other.append(compound)
        return(not_in_other)

    def _sanitise(self):
        """
        Remove all duplicate compounds and store the duplicates in the synonyms list
        """
        sanitised_list = []
        full_length = len(self.compound_list)
        for compound in self.compound_list:
            if not compound.in_list(sanitised_list, main_layer = self.main_layer):
                sanitised_list.append(compound)

        self.compound_list = sanitised_list
        self.sanitised = True

    def __repr__(self):
        """ Currently the same as str, allows for representation in terminal"""
        organism = self.organism
        delimiter = self.representation.delimiter
        rep = delimiter
        if organism is None:
            for compound in self.compound_list:
                rep = rep + str(compound) + delimiter
        else:
            chassis_state = organism
            for compound in self.compound_list:
                if chassis_state.compound_in_state(compound):
                    rep = rep + self.representation.color_begin + str(compound) + self.representation.printing_solved +  self.representation.color_end + delimiter
                else:
                    rep = rep + str(compound) + delimiter
        return(rep)

    def remove_cmpd_from_state(self, compound_to_remove):
        compound_inside_set = None
        for compound_index in range(len(self.compound_list)):
            if self.main_layer:
                if compound_to_remove.eq_main_layer(self.compound_list[compound_index]):
                    compound_inside_set = self.compound_list[compound_index]
                    break
            else:
                if compound_to_remove.eq_full_inchi_key(self.compound_list[compound_index]):
                    compound_inside_set = self.compound_list[compound_index]
                    break
        if not compound_inside_set is None:
            del self.compound_list[compound_index]  # delete the compound from the set

    def _naming(self):
        """Closely resembles str and repr"""
        organism = self.organism
        delimiter = '_'
        rep = ''
        for compound in self.compound_list:
            rep = rep + str(compound.name) + delimiter
        self.state_name = rep[:-1]

    def __str__(self):
        """ Currently the same as str, allows for representation in terminal"""
        organism = self.organism
        delimiter = self.representation.delimiter
        rep = delimiter
        if organism is None:
            for compound in self.compound_list:
                rep = rep + str(compound) + delimiter
        else:
            chassis_state = organism
            for compound in self.compound_list:
                if chassis_state.compound_in_state(compound):
                    rep = rep + self.representation.color_begin + str(compound) + self.representation.printing_solved +  self.representation.color_end + delimiter
                else:
                    rep = rep + str(compound) + delimiter
        return(rep)

    def clone(self):
        """
        Cloning the state is key for rollout.
        It allows for changing the sate without creating new nodes
        """
        return(ChemicalCompoundState(self.compound_list,
                                     organism = self.organism,
                                     representation = self.representation,
                                     available_rules = self.available_rules,
                                     main_layer = self.main_layer,
                                     biological_scorer = self.biological_scorer,
                                     chemical_scorer = self.chemical_scorer))

    def _check_sanitized_state(self):
        """
        Using the fact that state is ordered,
        sanitizing just verifies elements next to each other are different
        """
        sanitised = True
        try:
            if len(self.compound_list) >=2:
                for i in range(len(self.compound_list)-1):
                    if self.main_layer:
                        if self.compound_list[i].eq_main_layer(self.compound_list[i+1]):
                            # sanitised = False
                            raise ChemicalSanitisationException
                    else:
                        if self.compound_list[i].eq_full_inchi_key(self.compound_list[i+1]):
                            # sanitised = False
                            raise ChemicalSanitisationException
        except ChemicalSanitisationException as e:
            sanitised = False
            self.logger.error(e)
        return(sanitised)

    def compound_in_state(self, compound):
        """
        Verifies an element is in the set.
        """
        compound_in_state = False
        for i in range(0, len(self.compound_list)):
            if self.main_layer:
                compound_equal = (compound.eq_main_layer(self.compound_list[i]))
            else:
                compound_equal = (compound.eq_full_inchi_key(self.compound_list[i]))

            compound_in_state = compound_in_state or compound_equal
        return(compound_in_state)

    def _add_compound_to_state(self, compound):
        """
        Adding a compound to the state
        - if it is no already in it
        - using append then sorted according to https://groups.google.com/forum/#!topic/comp.lang.python/JZGNaCZnNR4
        - could be modified using bisect
        """
        if not self.compound_in_state(compound):
            self.compound_list.append(compound)
            self.compound_list = sorted(self.compound_list)

    def merge_states(self, state):
        """
        Adding a compound to the state
        - if it is no already in it
        - using append then sorted according to https://groups.google.com/forum/#!topic/comp.lang.python/JZGNaCZnNR4
        - could be modified using bisect
        """
        added = 0
        for compound in state.compound_list:
            if not self.compound_in_state(compound):
                added = added + 1
                self.compound_list.append(compound)
                self.compound_list = sorted(self.compound_list)
        self.logger.info("Added {} compounds when merging states".format(added))

    def ApplyMove(self, move, return_compounds = False):
        """
        Move should allow application with only the compound number and rule_ID.
        Then it can be applied through querrying and chemistry equally
        Move will be a Move object
        """
        compound_to_change = None
        for compound_index in range(len(self.compound_list)):
            if self.main_layer:
                if self.compound_list[compound_index].InChIKey.split("-")[0] == move.compound_id.split("-")[0]:
                    compound_to_change = self.compound_list[compound_index]
                    break
            else:
                if self.compound_list[compound_index].InChIKey == move.compound_id:
                    compound_to_change = self.compound_list[compound_index]
                    break
        product_set = compound_to_change.apply_transformation_with_move(move = move)
        if return_compounds:
            products = []
        del self.compound_list[compound_index]  # delete the compound from the set
        for product in product_set:
            # create compound:
            compound_to_add = product
            # add it to the state
            self._add_compound_to_state(compound_to_add)
            if return_compounds:
                products.append(compound_to_add)
        if return_compounds:
            return(products)

    def GetMoves(self, top_x = 5, chemical_score = True, biological_score = True, extension = False):
        """
        GetMoves should only return compound_index and rule_ID, as these will be used for DB query.
        It returns the ordered moves with their scores
        Extension considers the case when a tree is rerun with wider number of allowed children.
        If both biological and chemical scores, sort by product
        """

        ordering_moves = {}
        all_moves = {}
        for compound_index in range(len(self.compound_list)):
            compound = self.compound_list[compound_index]
            if self.organism.compound_in_state(compound):
                pass  # Do not transform solved compounds
            else:
                compound.set_max_moves(top_x)
                moves = compound.obtain_applicable_transformation_with_move(available_rules = self.available_rules,
                                                                            chemical_scorer = self.chemical_scorer,
                                                                            biological_scorer = self.biological_scorer,
                                                                            extension = extension)
                for move in moves:
                    if chemical_score and biological_score:
                        ordering_moves[move.name] = move.chemical_score * move.biological_score
                    elif chemical_score:
                        ordering_moves[move.name] = move.chemical_score
                    elif biological_score:
                        ordering_moves[move.name] = move.biological_score
                    else:
                        self.logger.warning("There is no way to order moves specified for state {}.format(self)")
                    all_moves[move.name] = move
        top_x_moves = []
        for element in sorted(ordering_moves, key = lambda x:ordering_moves[x], reverse = True)[:top_x]:
            top_x_moves.append(all_moves[element])
        return(top_x_moves)

    def GetRolloutMoves(self):
        """
        Returns moves for Rollout explicitely.
        No need for ordering, as it is done by the Rollout Policy in the MCTS Node object
        """

        all_moves = []
        for compound_index in range(len(self.compound_list)):
            compound = self.compound_list[compound_index]
            if self.organism.compound_in_state(compound):
                pass  # Do not transform solved compounds
            else:
                moves = compound.obtain_applicable_transformation_with_move(available_rules = self.available_rules,
                                                                            chemical_scorer = self.chemical_scorer,
                                                                            biological_scorer = self.biological_scorer,)
                for move in moves:
                    all_moves.append(move)
        return(all_moves)

    def GetResults_from_Compounds(self, rewarding = Basic_Rollout_Reward):
        """
        WARNING: deprecated as inchikeys are faster.
        Getting reward:
        - rewarding.penalty if nothing in sink
        - rewarding.full_state_reward if everything in sink
        - fraction otherwise
        Remarks:
        - This implementation is probably slower
        - This would typically have to be changed for biosensors
        - More complex functions can and should be implemented
        """
        # TODO: define separate reward class
        chassis_metabolites = self.organism
        compounds_in_chassis = 0
        for compound in self.compound_list:
            if chassis_metabolites.compound_in_state(compound):
                compounds_in_chassis = compounds_in_chassis + 1
        if compounds_in_chassis == 0:
            reward = rewarding.penalty
        elif compounds_in_chassis == self.length:
            reward = rewarding.full_state_reward
        else:
            reward = compounds_in_chassis/self.length
        return(reward)

    def GetResults_from_InChI_Keys(self, rewarding = Basic_Rollout_Reward, main_layer = True):
        """
        Getting reward:
        - rewarding.penalty if nothing in sink
        - rewarding.full_state_reward if everything in sink
        - fraction otherwise
        Remarks:
        - This implementation is probably faster
        - This would typically have to be changed for biosensors
        - More complex functions can and should be implemented
        """
        if main_layer:
            chassis_metabolites = [compound.InChIKey.split('-')[0] for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            compounds_in_chassis = 0
            for compound in self.compound_list:
                InChI_key = compound.InChIKey.split('-')[0]
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    compounds_in_chassis = compounds_in_chassis + 1
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
        else:
            chassis_metabolites = [compound.InChIKey for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            compounds_in_chassis = 0
            for compound in self.compound_list:
                InChI_key = compound.InChIKey
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    compounds_in_chassis = compounds_in_chassis + 1
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
        if compounds_in_chassis == 0:
            reward = rewarding.penalty
        elif compounds_in_chassis == self.length:
            reward = rewarding.full_state_reward
        else:
            reward = compounds_in_chassis/self.length
        return(reward)

    def GetSupplement_from_InChI_Keys(self, main_layer = True):
        """
        The aim is to return a compound that could be a supplement if only one compound in unsolved
        Otherwise, return None
        """
        not_in_chassis_compounds = []
        if main_layer:
            chassis_metabolites = [compound.InChIKey.split('-')[0] for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            compounds_in_chassis = 0
            for compound in self.compound_list:
                InChI_key = compound.InChIKey.split('-')[0]
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    compounds_in_chassis = compounds_in_chassis + 1
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
                else:
                    not_in_chassis_compounds.append(compound)
        else:
            chassis_metabolites = [compound.InChIKey for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            compounds_in_chassis = 0
            for compound in self.compound_list:
                InChI_key = compound.InChIKey
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    compounds_in_chassis = compounds_in_chassis + 1
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
                else:
                    not_in_chassis_compounds.append(compound)

        if compounds_in_chassis == 0:
            return(None)
        elif compounds_in_chassis == self.length:
            assert len(not_in_chassis_compounds) == 0
            return(None)
        elif len(not_in_chassis_compounds) >1:
            return(None)
        else:
            assert len(not_in_chassis_compounds) == 1
            return(not_in_chassis_compounds[0])

    def GetResultsForBiosensors(self, rewarding = Basic_Rollout_Reward, main_layer = True):
        found = False
        if main_layer:
            chassis_metabolites = [compound.InChIKey.split('-')[0] for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            for compound in self.compound_list:
                InChI_key = compound.InChIKey.split('-')[0]
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    found = True
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
                    break
        else:
            chassis_metabolites = [compound.InChIKey for compound in self.organism.compound_list]
            chassis_metabolites_names = [compound.name for compound in self.organism.compound_list]
            for compound in self.compound_list:
                InChI_key = compound.InChIKey
                insertion_point = bisect.bisect_left(chassis_metabolites, InChI_key)
                if (insertion_point < len(chassis_metabolites) and InChI_key == chassis_metabolites[insertion_point]):
                    found = True
                    compound.add_synonym_by_name(chassis_metabolites_names[insertion_point])
                    break
        if found:
            reward = rewarding.full_state_reward
        else:
            reward = rewarding.penalty
        return(reward)


def __cli():
    """Command line interface.
    Only used to run tests at the moment
    """

    help = "Simple compound for MCTS. The CLI applies --rsmarts to --csmiles."
    logging.warning("No CLI is available for chemical_compounds_state")

if __name__ == "__main__":
    from chemical_compounds_state import *
    logging_level = logging.DEBUG
    logging.basicConfig(stream=sys.stderr,
                        level=logging_level,
                        datefmt='%d/%m/%Y %H:%M:%S',
                        format='%(asctime)s -- %(levelname)s -- %(message)s')
    __cli()
