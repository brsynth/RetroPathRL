"""
Contains the compound class that allows for rule obtention, application and scoring.
It is the base element of the chemical compound state class that contains multiple compounds.
Mathilde Koch, INRA, 2019
"""

# General utility packages
import logging
import sys
import time
import pickle
import argparse

# # Memory tracking - for development
# from pympler import tracker, classtracker

from multiprocessing import Pool, TimeoutError, current_process
# from multiprocessing import Pool, Queue, Process, TimeoutError, Array, JoinableQueue
#
# Chemistry
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import AllChem
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError
from utilities.reactor.cli import worker_match, worker_fire, kill, RuleConversionError
# MCTS move wrapper
from move import Move
# MCTS scoring functions:
from biological_scoring import BiologicalFullScoringH, BiologicalFullScoringRetroNoH, BiologicalFullScoringFwdNoH, BiologicalFullScoringRetroH, BiologicalFullScoringFwdH, RandomBiologicalScorer

from chemical_scoring import RandomChemicalScorer, SubstrateChemicalScorer, SubandprodChemicalScorer, ConstantChemicalScorer

# Data importation
from rule_sets_similarity import get_rules_and_score, full_rules_forward_H, full_rules_retro_H, full_rules_forward_no_H, \
    full_rules_retro_no_H
# Configuration
from config import *

if use_toxicity:
    from compound_scoring import toxicity_scorer

# The recursion limit was raised for pickling of Trees.
sys.setrecursionlimit(50000)

def worker_standardisation(kwargs):
    """
    Deprecated. Used when implementing parallelisation.
    Apply standardisation on a compound.
    """
    ans = standardize_chemical(kwargs["rd_mol"],
                                add_hs=kwargs["add_Hs"],
                                rm_stereo=kwargs["rm_stereo"],
                                heavy = kwargs["heavy_standardisation"])
    return ans

def worker_standard_results(kwargs):
    """Standardises results from a rule application"""
    ans = standardize_results(kwargs["ans"], add_hs=kwargs["add_Hs"], rm_stereo=kwargs["rm_stereo"])
    return(ans)


def _ask_DB(DB_CACHE, DB_REPLACE, rule_id, substrate_id):
    """
    Queries the DB for results when available.

    @param DB_CACHE (bool): switch to use a database to cache fire results
    @param DB_REPLACE (bool): switch that indicates if results have to be overwritten
    @param rule_id (str): reaction rule ID
    @param substrate_id (str): substrate ID
    @returns ({'found': str, 'list_rdmols': [], 'list_stoechiometry': []}): reply dictionary
    """
    if DB_CACHE and not DB_REPLACE:
        # Look for result
        document_id = make_document_id(rule_id=rule_id, substrate_id=substrate_id)
        document = CACHE_MGR.find(document_id)
        # If results are in DB
        if document is not None:
            # Rebuild rdmols & return
            list_rdmols, list_stoechiometry = rdmols_from_document(document, build_from="inchi", add_hs=True)
            return {
                'found': True,
                'list_rdmols': list_rdmols,
                'list_stoechiometry': list_stoechiometry
            }
    # In all other cases
    return {
        'found': False,
        'list_rdmols': [],
        'list_stoechiometry': []
    }


def _moves_from_rdmols(original_compound, rdmols, move, main_layer,
                       chemical_scorer, clean_up, stereo, max_moves, fire_timeout,
                       chemical_scoring_configuration,
                       rule_id, list_stoechiometry):
    """
    Generates Move objects from rdmols (results from rule application)
    """
    cleaned_rdmols = []
    list_of_moves = []
    if clean_up:
        list_stoechiometry= []
    else:
        list_stoechiometry = list_stoechiometry
    if not len(rdmols):
        return cleaned_rdmols, list_of_moves, list_stoechiometry
    else:
        original_substrates_list = move.original_substrates_list
        original_products_list_list = move.original_products_list_list
        for set_number in range(len(rdmols)):
            move_inside_set = move.clone()
            list_rdmols = rdmols[set_number]
            stoechiometry = {}
            results_as_compound_objects = []
            # Looking whether rdmol objects repeat themselves.
            cleaned_set_list = []
            for mol_obj in list_rdmols:
                compound = Compound(rdkit_obj = mol_obj, stereo = stereo,
                                    max_moves = max_moves, fire_timeout = fire_timeout,
                                    chemical_scoring_configuration = chemical_scoring_configuration)
                # WARNING: forcing standardisation again due to stereo problems. Advise not to use stereo
                # Stoechiometry information
                if compound.InChIKey not in stoechiometry.keys():
                    stoechiometry[compound.InChIKey] = 1
                if clean_up:
                    if compound.in_list(results_as_compound_objects, main_layer = main_layer):
                        stoechiometry[compound.InChIKey] = stoechiometry[compound.InChIKey] + 1
                    else:
                        results_as_compound_objects.append(compound)
                        cleaned_set_list.append(mol_obj)
                else:
                    results_as_compound_objects.append(compound)
                    cleaned_set_list.append(mol_obj)
            if original_compound.in_list(results_as_compound_objects, main_layer = main_layer):
                original_compound.logger.info("Self generating rule for compound {} and rule {}, set {}. Move is discarded in {}".format(original_compound.name, rule_id, set_number, current_process().name))
            else:
                # Add the current set to the cleaned rdmol list, it passed all filters.
                move_inside_set.set_set_number(set_number)
                if clean_up:
                    move_inside_set.set_stoechiometry(stoechiometry)
                    list_stoechiometry.append(stoechiometry)
                else:
                    move_inside_set.set_stoechiometry(list_stoechiometry[set_number])

                move_inside_set.set_product_list(results_as_compound_objects)
                chemical_score = chemical_scorer.calculate(compound =original_compound, rule = move_inside_set.rid,
                                                        products = move_inside_set.full_product_list(),
                                                        original_substrates_list = original_substrates_list,
                                                        original_products_list_list = original_products_list_list)
                move_inside_set.set_chemical_score(chemical_score)
                move_inside_set.delete_intermediate_chemical_score()
                cleaned_rdmols.append(cleaned_set_list)
                list_of_moves.append(move_inside_set)
        assert len(list_of_moves) == len(list_stoechiometry)
        return(cleaned_rdmols, list_of_moves, list_stoechiometry)

def calculate_from_available_rules(arguments):
    """
    Deprecated. Used with __obtain_move_parallel method of compound.
    """
    rule_id, rule_characteristics, substrate_info = arguments
    rdmol = substrate_info["rdmol"]
    inchikey = substrate_info["inchikey"]
    substrate_id = substrate_info["inchikey"]
    stereo = substrate_info["stereo"]
    max_moves = substrate_info["max_moves"]
    fire_timeout = substrate_info["fire_timeout"]
    cmp = substrate_info["cmp"]
    main_layer = substrate_info["main_layer"]
    chemical_scorer = substrate_info["chemical_scorer"]

    db_reply = _ask_DB(DB_CACHE=DB_CACHE, DB_REPLACE=DB_REPLACE, rule_id=rule_id, substrate_id=substrate_id)
    rdmols_from_DB = db_reply['list_rdmols']
    list_stoechiometry = db_reply['list_stoechiometry']

    # rdmols_from_DB = None
    if not db_reply['found']:
        # Initilise the rule
        rsmarts = rule_characteristics["Rule_SMARTS"]
        try:
            rd_rule = AllChem.ReactionFromSmarts(rsmarts)
            rd_rule.Initialize()
        except Exception as e:
            raise RuleConversionError(e) from e
        r = RuleBurnerCore(rd_rule = rd_rule, rd_mol = rdmol)
        ans = r.fire()
        if ans != ():
        # Creating the basic Move object for metadata
            move = Move(rsmart = rule_characteristics["Rule_SMARTS"],
                        rid = rule_id,
                        compound_id = inchikey,
                        biological_score = rule_characteristics["biological_score"],
                        EC_number= rule_characteristics["EC_number"])
            rdmols, failed = standardize_results(ans, add_hs=add_Hs, rm_stereo=not stereo)  # !!!! add_hs=False To be used only to fill the DB
            clean_up = True
        else:
            rdmols = None
    else:
        rdmols = rdmols_from_DB
        clean_up = False
        # Stoechiometry needs to come from the DB
        move = Move(rsmart = rule_characteristics["Rule_SMARTS"],
                    rid = rule_id,
                    compound_id = inchikey,
                    biological_score = rule_characteristics["biological_score"],
                    EC_number= rule_characteristics["EC_number"])

        # print("rdmols for rule {} is {}".format(rule_id, rdmols))
    if not rdmols is None:
        original_substrates_list = rule_characteristics["substrate_ECFP"]
        original_products_list_list = rule_characteristics["products_ECFP"]
        move.set_intermediate_chemical_score(original_substrates_list, original_products_list_list)
        cleaned_rdmols, list_of_moves, list_stoechiometry = _moves_from_rdmols(original_compound = cmp,
                                                            rdmols = rdmols,
                                                            move = move,
                                                            main_layer =main_layer,
                                                            chemical_scorer = chemical_scorer,
                                                            clean_up = clean_up,
                                                            stereo = stereo,
                                                            max_moves = max_moves,
                                                            fire_timeout = fire_timeout,
                                                            rule_id = rule_id,
                                                            list_stoechiometry = list_stoechiometry)
        logging.info("There are list_of_moves {} for rule {}".format(list_of_moves, rule_id))
    else:
        cleaned_rdmols = None
        list_of_moves = []
    if DB_CACHE and DB_REPLACE:
        # Format document
        # Here, check it's not empty.
        if not cleaned_rdmols is None:
            list_list_inchikeys, list_list_inchis, list_list_smiles = handle_results(cleaned_rdmols)
        else:
            list_list_inchikeys, list_list_inchis, list_list_smiles, list_stoechiometry = [], [], [], []
        document = as_document(
                rule_id=rule_id,
                substrate_id=substrate_id,
                list_list_inchikeys=list_list_inchikeys,
                list_list_inchis=list_list_inchis,
                list_list_smiles=list_list_smiles,
                list_stoechiometry = list_stoechiometry)
        # print("Exported document to dB {}".format(document))
        # Insert
        CACHE_MGR.insert(document)
    # logging.info("For rule {}, return moves {}".format(rule_id, list_of_moves))
    # print("For rule {}, return moves {}".format(rule_id, list_of_moves))
    return(list_of_moves)
    # time.sleep(random.randint(1, 10))
    # return(rule_id)

def unpickle(file_name, type = "compound", folder_address = "pickled_data"):
    with open('{}/{}_{}.pkl'.format(folder_address, type, file_name), 'rb') as input:
        return(pickle.load(input))

def worker_cmp(worker_id, in_queue, out_queue, factory_task, timeout=1):
    """
    DEPRECATED
    One worker handle one renewable pool running one process"""
    if DB_CACHE:
        logging.info("Using DB cache for worker {}".format(worker_id))
        CACHE_MGR = Manager(replace=DB_REPLACE)
        CACHE_MGR.connect()
    pool_version = 0
    pool_worker = None
    while True:
        if pool_worker is None:
            pool_worker = Pool(processes=1)
            pool_version = pool_version + 1
        # print(killing_array[int(worker_id)])
        if killing_array[int(worker_id)] > 0:
            logging.debug("Flushing the queue with worker {}".format(worker_id))
            while not in_queue.empty():
                in_queue.get()
                # in_queue.task_done()
            kill(pool_worker)
            break
            return()

        in_value = in_queue.get()
        # logging.debug("ID: {}, with rule ID {}".format(worker_id, in_value[0]))
        if in_value == POISON_PILL:
            print("{}'s speaking: Oh nooo nooo!!".format(worker_id))
            kill(pool_worker)
            break
        # Do something smarts
        ares = pool_worker.apply_async(factory_task, (in_value,))
        try:
            ans = ares.get(timeout=timeout)
            if ans!= []:  #TODO adapt it for other types of respones from factories
                out_queue.put(ans)
            # out_queue.put(ans)

            # print(ans)
        except TimeoutError as e:
            kill(pool_worker)
            pool_worker = None
            logging.info("Worker {}. Pool {} has been killed for Timeout".format(worker_id, pool_version))
    return()


class FinishedException(Exception):
    def __init__(self, msg = "Reached allowed number of transformations"):
        self._msg = msg

    def __str__(self):
        return self._msg

class CompoundException(Exception):
    """Home made exception for Compound Class."""
    pass

class UnSpecifiedTransformationError(CompoundException):
    """Raised when a Transformation is not specified by either smart or RuleID."""

    def __init__(self, msg = "UnSpecifiedTransformationError"):
        self._msg = msg

    def __str__(self):
        return self._msg

class CompoundDefinitionException(CompoundException):
    """When compound structure is not properly defined"""
    def __init__(self, msg = "Not properly specified, lacking smiles or InChI"):
        self._msg = msg

    def __str__(self):
        return self._msg

class Compound(object):
    """
    Basic chemical compound object for MCTS.
    Base element of the chemical compound state.
    This is where
    - the chemistry is performed (obtain reaction rules that apply)
    - the scoring is performed
    The higher function objects such as Tree and MCTS Node depend on this object class

    """
    logger = logging.getLogger(__name__)

    def __init__(self,
                 csmiles = None,
                 InChI = None,
                 name = None,
                 stereo = False,
                 rdkit_obj = None,
                 synonyms = None,
                 synonyms_names = None,
                 max_moves = 5,
                 single_worker = True,
                 fire_timeout = 1,
                 standardisation_timeout = 5,
                 heavy_standardisation = False,
                 force_add_H = None,
                 chemical_scoring_configuration = {"biological_score_cut_off": 0.1,
                 "substrate_only_score_cut_off": 0.3,
                 "chemical_score_cut_off": 0.3},
                 toxicity = None):
        """
        Initialising a chemical coumpound.
        Mostly an RDKit wrapper with
            - biological score
            - chemical score
        The name is only for the visualisation functions - set to inchikey if not specficied
        """
        self.stereo = stereo
        self.force_add_H = force_add_H
        self.single_worker = single_worker
        if self.single_worker:
            self._pool = None
        if rdkit_obj is None:
            if InChI is None and csmiles is None:
                raise CompoundDefinitionException
            else:
                self.InChI = InChI
                self.smiles = csmiles
                self._rdmol_set_up(heavy_standardisation)  # Standardise
            self._calculate_InChIKey()
            if name is None:
                self.name = self.InChIKey
            else:
                self.name = name
            if heavy_standardisation:
                self.logger.debug("Used heavy standardisation for {}".format(self.name))
        else:
            self.rdmol = rdkit_obj
            self.csmiles = Chem.MolToSmiles(self.rdmol,isomericSmiles=True, allHsExplicit=False)
            self._calculate_InChIKey()

            if heavy_standardisation:
                self.logger.debug("Using heavy standardisation for {}".format(self.name))
            if name is None:
                self.name = self.InChIKey
            else:
                self.name = name
            # self.logger.debug("Considering mol file to be standardised for {}".format(self))
        self._rdkit_logger = RDLogger.logger()
        if synonyms is None:
            self.synonyms = [self.csmiles]
        else:
            self.synonyms = synonyms
        if synonyms_names is None:
            self.synonyms_names = [self.name]
            if self.name != self.InChIKey:
                self.synonyms_names.append(self.InChIKey)
        else:
            self.synonyms_names = synonyms_names
        self._fire_timeout = fire_timeout
        self.chemical_scoring_configuration = chemical_scoring_configuration
        self.standardisation_timeout = standardisation_timeout
        self.max_moves = max_moves
        self.found_moves = None
        self.ECFP = None
        self.toxicity = toxicity
        if use_toxicity:
            if self.toxicity is None:
                self.toxicity = min(toxicity_scorer.calculate(self), 0)
                self.logger.debug("Compound {} has a toxicity score of {}".format(self, self.toxicity))


    def eq_full_inchi_key(self, other):
        """
        Two compounds are identical if their full InChI keys are identical
        Remarks:
        - necessitates 'slow' InChIKey calculation
        - allows for fast comparison between compounds
        - only works for exactly macthing compounds
        """
        cs_InChIkey_self = self.InChIKey
        cs_InChIkey_other = other.InChIKey
        return cs_InChIkey_self == cs_InChIkey_other  # Assuming no InChIKey collisions

    def eq_main_layer(self, other):
        """
        Two compounds are identical if their InChI keys are identical
        Remarks:
        - necessitates 'slow' InChIKey calculation
        - allows for faster comparison between compounds
        - works for compounds that have rather different standardisation (necessary in RP2 and here)
        - will be used mostly for chassis, can also be used for speed up
        """
        cs_InChIkey_self = self.InChIKey.split("-")[0]
        cs_InChIkey_other = other.InChIKey.split("-")[0]
        return cs_InChIkey_self == cs_InChIkey_other  # Assuming no InChIKey collisions

    def __hash__(self):
        """
        Hash method should return an integer
        """
        return(hash(self.InChIKey))

    def in_list(self, list_of_compounds, main_layer = True):
        in_list = False
        for compound_in_list in list_of_compounds:
            if main_layer:
                equality = self.eq_main_layer(compound_in_list)
                if equality:
                    in_list = True
                    compound_in_list.add_synonym(self)
                    break
            else:
                equality = self.eq_full_inchi_key(compound_in_list)
                if equality:
                    in_list = True
                    compound_in_list.add_synonym(self)
                    break
        return(in_list)

    def __lt__(self, other):
        """
        Defines an order on compounds
        Used for faster look-up in ordered lists to see if a compound belongs to a set
        Current implementation: alphabetical order on InChIKeys
        """
        cs_InChIkey_self = self.InChIKey
        cs_InChIkey_other = other.InChIKey
        return (cs_InChIkey_self < cs_InChIkey_other)

    def __str__(self):
        """Compound representation is its name"""
        return(str(self.name))

    def __repr__(self):
        """Compound representation is its name"""
        return(str(self.name))

    def detailed_print(self):
        """Compound detailed information is printed"""
        text_bloc_start = 'Compound {} informations:\n'.format(self.name)
        synyms = "Its synonyms are {} \n".format(self.synonyms_names)
        inchi = "Its InChI is {}\n".format(self.InChI)
        smile = "Its canonical smiles is {}\n".format(self.csmiles)
        return(text_bloc_start + synyms + inchi + smile)

    def set_max_moves(self, max_moves):
        """Notable use include tree extension"""
        self.max_moves = max_moves

    def set_chemical_scoring_configuration(self, chemical_scoring_configuration):
        """Notable use include tree extension"""
        self.chemical_scoring_configuration = chemical_scoring_configuration

    def save(self, file_name = None, folder_address = "pickled_data"):
        if file_name is None:
            file_name = self.name
        file_saving = open('{}/compound_{}.pkl'.format(folder_address, file_name), 'wb')
        pickle.dump(self, file_saving)

    def clone(self):
        """
        Cloning is tested in the test files.
        Will need to include moves when they become an attribute. That should become much faster.
        """
        new_smiles = self.csmiles
        try:
            toxicity = self.toxicity
        except AttributeError:
            toxicity = None
        return(Compound(rdkit_obj = self.rdmol,
                        name = self.name,
                        stereo = self.stereo,
                        synonyms = self.synonyms,
                        synonyms_names =self.synonyms_names,
                        max_moves = self.max_moves,
                        single_worker = self.single_worker,
                        fire_timeout = self._fire_timeout,
                        chemical_scoring_configuration = self.chemical_scoring_configuration,
                        standardisation_timeout = self.standardisation_timeout,
                        toxicity = toxicity))

    def naming(self, name):
        """Allows for renaming, mostly when a state is initialised with a list"""
        self.name = name

    def add_synonym(self, compound):
        """
        Adds a synonym to this compound.
        When another compound was deemed equal to compound.
        - Currently saving it as a smiles string to be able find the structure
        """
        if compound.csmiles not in self.synonyms:
            self.synonyms.append(compound.csmiles)
            self.synonyms_names.append(compound.name)
    def add_synonym_by_name(self, new_name):
        """
        Adds a synonym to this compound.
        When another compound was deemed equal to compound.
        - Currently saving it as a smiles string to be able find the structure
        """
        if new_name not in self.synonyms_names:
            self.synonyms_names.append(new_name)

    def _rdmol_set_up(self, heavy_standardisation):
        """
        Set up the RDKit underlying object. Taken from Thomas.
        Allows for setting up of the RDKit object that allows for
        - ECFP calculation if need be
        - InChIKey calculation (necessary with current implementation)
        CANNOT BE SIMPLER THAN THAT FOR OTHER CODE MPORTATION
        """
        if self.smiles is not None:
            try:
                rd_mol = Chem.MolFromSmiles(self.smiles, sanitize=False)  # Important: Sanitize = False
                rd_mol = self._standardize_chemical(rd_mol, heavy_standardisation)
            except Exception as e:
                raise ChemConversionError(e) from e
        else:
            try:
                rd_mol = Chem.inchi.MolFromInchi(self.InChI, sanitize=False)  # Important: Sanitize = False
                rd_mol = self._standardize_chemical(rd_mol, heavy_standardisation)
            except Exception as e:
                raise ChemConversionError(e) from e
        self.rdmol = rd_mol
        self.csmiles = Chem.MolToSmiles(self.rdmol,isomericSmiles=True, allHsExplicit=False)

    def _calculate_InChIKey(self):
        """
        Remark: we are not generating standard InChI keys or InChIs.
        It's okay as long as we're consistent
        """
        self.InChI = Chem.rdinchi.MolToInchi(self.rdmol)[0]
        self.InChIKey = Chem.rdinchi.MolToInchiKey(self.rdmol)

    def _standardize_chemical(self, rd_mol, heavy_standardisation = False):
        """
        Standardize a chemical using RDKit sanitize method.
        Taken from rule_fire RuleBurner class

        :param      rd_mol:  RDKit mol object
        :returns    nothing, in place operation
        """
        rm_stereo = not self.stereo
        if self.force_add_H is None:
            return(standardize_chemical(rd_mol, add_hs=add_Hs, rm_stereo=rm_stereo, heavy = heavy_standardisation))
        else:
            return(standardize_chemical(rd_mol, add_hs=self.force_add_H, rm_stereo=rm_stereo, heavy = heavy_standardisation))

    def _run_with_timeout(self, worker, kwargs, timeout=1):
        """Generic wrapper making use of multiprocessing to garantee effective timeout.

        :param  worker:     function to be called
        :param  kwargs:     dictionnary of args to be passed to the called function
        :param  timeout:    int, timeout
        :returns            depends of the worker function
        taken from cli of reactor.
        I might want to change this for my usage if I do more parallel things than this.
        """
        if self._pool is None:
            self._pool = Pool(processes=4)
        try:
            start_time = time.time()
            res = self._pool.apply_async(worker, kwds={'kwargs': kwargs})
            ans = res.get(timeout=timeout)
            end_time = time.time()
            exec_time = round(end_time - start_time, 4)
        except TimeoutError as e:
            kill(self._pool)
            self._pool = Pool(processes=1)
            raise e
        except Exception as e:
            raise e

        return ans, exec_time

    def _get_Biological_scorer(self, biological_scorer):
        if biological_scorer == "RandomBiologicalScorer":
            policy = RandomBiologicalScorer
        elif biological_scorer == "BiologicalFullScoring":
            policy = BiologicalFullScoring
        elif biological_scorer == "BiologicalFullScoringsim":
            policy = BiologicalFullScoringsim
        elif biological_scorer == "BiologicalBiosensorScoring":
            policy = BiologicalBiosensorScoring
        elif biological_scorer == "BiologicalFullScoringRetroH":
            policy = BiologicalFullScoringRetroH
        elif biological_scorer == "BiologicalFullScoringFwdH":
            policy = BiologicalFullScoringFwdH
        elif biological_scorer == "BiologicalFullScoringRetroNoH":
            policy = BiologicalFullScoringRetroNoH
        elif biological_scorer == "BiologicalFullScoringFwdNoH":
            policy = BiologicalFullScoringFwdNoH
        elif biological_scorer == "BiologicalFullScoringH":
            policy = BiologicalFullScoringH
        else:
            raise NotImplementedError(biological_scorer)
        return(policy)

    def _get_Chemical_scorer(self, chemical_scorer):
        if chemical_scorer == "RandomChemicalScorer":
            policy = RandomChemicalScorer
        elif chemical_scorer == "SubstrateChemicalScorer":
            policy = SubstrateChemicalScorer
        elif chemical_scorer == "SubandprodChemicalScorer":
            policy = SubandprodChemicalScorer
        elif chemical_scorer == "ConstantChemicalScorer":
            policy = ConstantChemicalScorer
        else:
            raise NotImplementedError(chemical_scorer)
        return(policy)

    def _obtain_available_rules(self, available_rules, chemical_scorer, biological_scorer):
        """
        The aim of this function is to obtain all rules that can potentially be applied.
        It has options (set in the chemical_scoring_configuration attribute) for:
        - using a biological score cut off
        - using a substrate noly cut-off
        - Warning: full chemical score cut-off would make no sense here, it will be after rule application
        """
        rules_with_cut_off = {}
        # Extracting rules with a threshold
        if not self.found_moves is None:
            already_tested_ids = [move.rid for move in self.found_moves]
        else:
            already_tested_ids = []

        for rule_id, rule_characteristics in available_rules.items():
            if rule_id not in already_tested_ids:
                if float(rule_characteristics['biological_score']) >= self.chemical_scoring_configuration['biological_score_cut_off']:
                    try:
                        original_substrates_list = rule_characteristics["substrate_ECFP"]
                        if chemical_scorer is SubstrateChemicalScorer or chemical_scorer is SubandprodChemicalScorer:
                            substrate_only_score = SubstrateChemicalScorer.calculate(self, products = None, rule = rule_id, original_substrates_list = original_substrates_list, original_products_list_list = None)
                        else:
                            substrate_only_score = 1
                    except KeyError:
                        substrate_only_score = 1
                        rule_characteristics["substrate_ECFP"] = None
                        rule_characteristics["products_ECFP"] = None
                    if substrate_only_score >= self.chemical_scoring_configuration["substrate_only_score_cut_off"]:
                        rules_with_cut_off[rule_id] = rule_characteristics
                        rules_with_cut_off[rule_id]["substrate_chemical_score"] = substrate_only_score

        # Sorting rules
        sorted_rules = sorted(rules_with_cut_off.items(), key=lambda item: (item[1]["biological_score"] * item[1]["substrate_chemical_score"], item[0]), reverse=True)
        return(sorted_rules)

    def _apply_rule(self, rule_id, rd_rule):
        application_time_start = time.time()
        self.logger.debug("Applying rule {} on {} -{} at {}".format(rule_id, self, self.csmiles, time.time()))
        if self.single_worker:
            try:
                kwargs = {
                        'rd_rule': rd_rule,
                        'rd_mol':  self.rdmol
                        }
                ans, fire_exec_time = self._run_with_timeout(
                        worker=worker_fire, kwargs=kwargs,
                        timeout=self._fire_timeout
                        )
            except ChemConversionError as e:
                ans = None
                fire_error = str(e)
                logging.warning(e)
            except TimeoutError as e:
                ans = None
                fire_error = str(e)
                logging.error('TIMEOUT: cid={}, rid={}'.format(self, rule_id))
                # time.sleep(0.01)
            except Exception as e:
                ans = None
                fire_error = str(e)
                logging.warning(e)
        else:
            r = RuleBurnerCore(rd_rule = rd_rule, rd_mol = self.rdmol)
            ans = r.fire()
        round_time = time.time() - application_time_start

        self.logger.debug("Applied rule {} on {} -{} in {}".format(rule_id, self, self.csmiles, round_time))
        return(ans)

    def obtain_applicable_transformation_with_move(self,
                                                   available_rules,
                                                   chemical_scorer = "RandomChemicalScorer",
                                                   biological_scorer = "RandomBiologicalScorer",
                                                   main_layer = True,
                                                   extension = False):
        """
        Used to find applicable rules that respect scoring configuration.
        Calls functions for
        - obtaining rules
        - applying individual rules (or SQL)
        - scoring
        """
        #
        start_time = time.time()
        biological_scorer = self._get_Biological_scorer(biological_scorer)
        chemical_scorer = self._get_Chemical_scorer(chemical_scorer)
        self.logger.info("Starting move obtention process with compound {} ({}) at {}".format(self, self.csmiles, start_time))

        if use_cache and self.InChIKey in home_made_cache.keys() and not extension:
            # If already calculated during this run of the tree, and not extenting moves.
            matching_moves = home_made_cache[self.InChIKey]
            self.logger.info("Retrieved moves in {} for compound {} - with caching".format(time.time()- start_time, self.InChIKey))
        else:
            matching_moves = []
            # Obtain rules sorted by the score I'm interested in.
            sorted_rules = self._obtain_available_rules(available_rules, chemical_scorer, biological_scorer)
            self.logger.info("{} available rules".format(len(sorted_rules)))
            for rule_id, rule_characteristics in sorted_rules:
                # Those rules are sorted by biological * chemical susbtrate score, maybe with cut offs.
                # If reach number of moves, stop looking
                if len(matching_moves) >= self.max_moves:
                    break

                # Give a chance with the DB cache
                db_reply = _ask_DB(DB_CACHE=DB_CACHE, DB_REPLACE=DB_REPLACE, rule_id=rule_id, substrate_id=self.InChIKey)

                if not db_reply['found']:
                    # Force calculation of the move with rdkit and chemistry
                    # Initilise the rdkit rule
                    start_application_time = time.time()
                    rsmarts = rule_characteristics["Rule_SMARTS"]
                    try:
                        rd_rule = AllChem.ReactionFromSmarts(rsmarts)
                        rd_rule.Initialize()
                    except Exception as e:
                        raise RuleConversionError(e) from e
                    # Apply the rule
                    ans = self._apply_rule(rule_id, rd_rule)
                    if not ans is None:
                        move = Move(rsmart = rule_characteristics["Rule_SMARTS"],
                                        rid = rule_id,
                                        compound_id = self.InChIKey,
                                        biological_score = rule_characteristics["biological_score"],
                                        EC_number= rule_characteristics["EC_number"])
                        kwargs_standard = {
                              "ans": ans,
                              "add_Hs": add_Hs,
                              "rm_stereo": not self.stereo
                              }
                        try:
                            ans, exec_time = self._run_with_timeout(worker_standard_results, kwargs_standard, timeout=self.standardisation_timeout)
                            rdmols, failed = ans
                        except TimeoutError as e:
                            self.logger.warning(e)
                            rdmols = []
                        #  rdmols, failed = standardize_results(ans, add_hs=add_Hs, rm_stereo=not self.stereo)  # !!!! add_hs=False To be used only to fill the DB
                        clean_up = True
                        list_stoechiometry = []
                    else:
                        rdmols = []
                        list_stoechiometry = []
                else:
                    rdmols = db_reply['list_rdmols']
                    list_stoechiometry = db_reply['list_stoechiometry']
                    clean_up = False
                    move = Move(rsmart = rule_characteristics["Rule_SMARTS"],
                                rid = rule_id,
                                compound_id = self.InChIKey,
                                biological_score = rule_characteristics["biological_score"],
                                EC_number= rule_characteristics["EC_number"])
                # If rdmols is still empty, then it's because the rule doesn't apply
                original_substrates_list = rule_characteristics["substrate_ECFP"]
                original_products_list_list = rule_characteristics["products_ECFP"]
                try:
                    move.set_intermediate_chemical_score(original_substrates_list, original_products_list_list)
                    cleaned_rdmols, list_of_moves, list_stoechiometry = _moves_from_rdmols(original_compound = self,
                                                                        rdmols = rdmols,
                                                                        move = move,
                                                                        main_layer= main_layer,
                                                                        chemical_scorer = chemical_scorer,
                                                                        clean_up = clean_up,
                                                                        stereo = self.stereo,
                                                                        max_moves = self.max_moves,
                                                                        fire_timeout = self._fire_timeout,
                                                                        chemical_scoring_configuration = self.chemical_scoring_configuration,
                                                                        rule_id = rule_id,
                                                                        list_stoechiometry = list_stoechiometry)
                    for move_inside_set in list_of_moves:
                        if not move_inside_set.in_list(list_moves = matching_moves, main_layer= main_layer):
                            # Add the move cut off on chemical score here, as I'll still want to have the rule in the DB.
                            if move_inside_set.chemical_score >= self.chemical_scoring_configuration["chemical_score_cut_off"]:
                                move_inside_set.calculate_rsmiles(self)
                                matching_moves.append(move_inside_set)
                    end_time = time.time()

                    if DB_CACHE and (DB_REPLACE or not db_reply['found']) and (end_time - start_time >=DB_time): # could add something about time for execution here
                        # Format document
                        # Here, check it's not empty.
                        list_list_inchikeys = [move.product_list for move in list_of_moves]
                        list_list_inchikeys, list_list_inchis, list_list_smiles = handle_results(cleaned_rdmols)
                        document = as_document(
                                rule_id=rule_id, substrate_id=self.InChIKey,
                                list_list_inchikeys=list_list_inchikeys,
                                list_list_inchis=list_list_inchis,
                                list_list_smiles=list_list_smiles,
                                list_stoechiometry = list_stoechiometry)
                        # print("Exported document to dB {}".format(document))
                        # Insert
                        CACHE_MGR.insert(document)
                        # total_clean_up = total_clean_up + time.time() - clean_up_start
                except UnboundLocalError as e:
                    self.logger.critical("Error {} when trying rule {} on compound {}. Moving to next rule".format(e, rule_id, self), exc_info=True)
            self.logger.info("For compound {}, matching moves are {}".format(self, matching_moves))
            if use_cache:
                self.logger.info("Caching results for compound {}. The cache contained {}".format(self.InChIKey, len(home_made_cache.keys())))
                home_made_cache[self.InChIKey] = matching_moves
            self.logger.info("Calculated moves in {} for compound {} - without caching".format(time.time()- start_time, self.InChIKey))

        # self.logger.debug("For compound {}, time spent in applying is: {}, sandardisation: {}, cleanup: {}".format(self, total_rule_application, total_standardisation, total_clean_up))
        if self.single_worker and not self._pool is None:
            # Try this to see if I still have a bug
            self._pool.close()
            self._pool = None
            # 
            # Should have done a cloe, but for pickling doing None.
        self.found_moves = matching_moves
        return(matching_moves)

    def __obtain_move_parallel(self, available_rules,
                            chemical_scorer = "RandomChemicalScorer",
                            biological_scorer = "RandomBiologicalScorer",
                            main_layer = True,
                            extension = False,
                            num_workers = 4,
                            timeout = 5,
                            individual_cmp_time = 120):
        """
        DEPRECATED.
        Used to develop parallelisation.
        Current issues involve the filling/emptying of the queue.
        Used to find applicable rules.
        Calls functions for
        - obtaining rules
        - applying individual rules (or DB_CACHE)
        - scoring
        """
        #
        start_time = time.time()
        biological_scorer = self._get_Biological_scorer(biological_scorer)
        chemical_scorer = self._get_Chemical_scorer(chemical_scorer)
        self.logger.info("Starting move obtention process with compound {} ({}) at {}".format(self, self.csmiles, start_time))
        self.logger.info("Defining parallel tools")
        in_queue = Queue()
        out_queue = Queue()
        global killing_array
        killing_array = Array("d", num_workers)
        factory_worker = calculate_from_available_rules
        global POISON_PILL
        POISON_PILL = "STOP"
        # Setting up the team
        team = list()
        for i in range(0, num_workers):
            team.append(Process(target=worker, kwargs={'worker_id': i,'in_queue': in_queue,
                                                        'out_queue': out_queue, 'timeout': timeout,
                                                        "factory_task": factory_worker,
                                                        "killing_array": killing_array, "POISON_PILL": POISON_PILL}))

        if use_cache and self.InChIKey in home_made_cache.keys() and not extension:
            # If already calculated during this run of the tree, and not extenting moves.
            matching_moves = home_made_cache[self.InChIKey]
            self.logger.info("Retrieved moves in {} for compound {} - with caching".format(time.time()- start_time, self.InChIKey))
        else:
            sorted_rules = self._obtain_available_rules(available_rules)
            substrate_info = {"cmp": self,
                            "rdmol": self.rdmol,
                            "inchikey": self.InChIKey,
                            "stereo": self.stereo,
                            "max_moves": self.max_moves,
                            "fire_timeout": self._fire_timeout,
                            "main_layer": main_layer,
                            "chemical_scorer": chemical_scorer}
            task_list = [(rule_id, rule_characteristics, substrate_info) for rule_id, rule_characteristics in sorted_rules]
            self.logger.info("Adding elements to the queue")
            #
            for task in task_list:
                in_queue.put(task)
            self.logger.info("{} available rules".format(len(sorted_rules)))

            self.logger.info("Starting the team")
            for team_player in team:
                team_player.start()
            matching_moves = []
            waiting_loops = 0
            try:
                while True:
                    try:
                        if len(matching_moves) >= self.max_moves:
                            raise FinishedException
                        list_of_moves = out_queue.get()
                        if list_of_moves != []:
                            #     logging.info("Returned something empty {}".format(list_of_moves))
                            # else:
                            for move_inside_set in list_of_moves:
                                if not move_inside_set.in_list(list_moves = matching_moves, main_layer= main_layer):
                                    matching_moves.append(move_inside_set)
                            # Will need to check here whether it's good or not
                            logging.info("total answer currently is {}".format(len(matching_moves)))
                    except Empty:
                        # time.sleep(timeout)
                        # print("nothing new")
                        waiting_loops = waiting_loops + 1

                        # print('Currently empty')
                    # print(len(total_answers))
                    finished_moves = len(matching_moves) >= self.max_moves
                    finished_input_queue = in_queue.empty()
                    finished_output_queue = out_queue.empty()
                    if finished_input_queue or finished_output_queue:
                        time.sleep(timeout)
                        print("Extending time because finished input {} or finished output {}".format(finished_input_queue, finished_output_queue))
                    finished_input_queue = in_queue.empty()
                    finished_output_queue = out_queue.empty()
                    if finished_moves or (finished_input_queue and finished_output_queue) or (time.time()- start_time > individual_cmp_time):
                        print("stopped because: lengh {}, input {}, output {}, time {}".format(finished_moves, finished_input_queue, finished_output_queue, time.time()- start_time > individual_cmp_time))
                        print("Found {} moves".format(len(matching_moves)))
                        raise FinishedException
                        # break
            except FinishedException:
                self.logger.info("Finished the number of moves")
                self.logger.info("Calculated moves in {} for compound {} - without caching".format(time.time()- start_time, self.InChIKey))

                # print("using killing array")

                for i in range(len(team)):
                    killing_array[i] = 1
                in_queue = Queue()
                print("other killing streatgy")
                # in_queue.close()
                # in_queue.join_thread()
                while any([team_player.is_alive() for team_player in team]):
                    in_queue.put(POISON_PILL)
                in_queue.join()
                # time.sleep(0.01)
                # for team_player in team:
                #     team_player.join()
                #     print("joining")
                for team_player in team:
                    if team_player.is_alive():
                        self.logger.warning("Some workers are still alive")

        self.logger.info("For compound {}, matching moves are {}".format(self, matching_moves))
        if use_cache:
            self.logger.info("Caching results for compound {}. The cache contained {}".format(self.InChIKey, len(home_made_cache.keys())))
            home_made_cache[self.InChIKey] = matching_moves
        self.logger.info("Calculated moves in {} for compound {} - without caching".format(time.time()- start_time, self.InChIKey))

        # print(matching_moves)
        return(matching_moves)

    def __obtain_move_parallel_with_factory(self, available_rules,
                            chemical_scorer = "RandomChemicalScorer",
                            biological_scorer = "RandomBiologicalScorer",
                            main_layer = True,
                            extension = False,
                            num_workers = 4,
                            timeout = 5,
                            individual_cmp_time = 120,
                            tasks_per_chunk = 1000):
        """
        DEPRECATED. Used for parallelisation.
        Used to find applicable rules.
        Calls functions for
        - obtaining rules
        - applying individual rules (or SQL)
        - scoring
        """
        #
        start_time = time.time()
        biological_scorer = self._get_Biological_scorer(biological_scorer)
        chemical_scorer = self._get_Chemical_scorer(chemical_scorer)
        self.logger.info("Starting move obtention process with compound {} ({}) at {}".format(self, self.csmiles, start_time))

        if use_cache and self.InChIKey in home_made_cache.keys() and not extension:
            # If already calculated during this run of the tree, and not extenting moves.
            matching_moves = home_made_cache[self.InChIKey]
            self.logger.info("Retrieved moves in {} for compound {} - with caching".format(time.time()- start_time, self.InChIKey))
        else:
            sorted_rules = self._obtain_available_rules(available_rules)
            substrate_info = {"cmp": self,
                            "rdmol": self.rdmol,
                            "inchikey": self.InChIKey,
                            "stereo": self.stereo,
                            "max_moves": self.max_moves,
                            "main_layer": main_layer,
                            "chemical_scorer": chemical_scorer}
            task_list = [(rule_id, rule_characteristics, substrate_info) for rule_id, rule_characteristics in sorted_rules]
            self.logger.info("Adding elements to the queue")
            #
            chunks = round(len(task_list)/tasks_per_chunk) + 1
            task_chunks = {}
            indices_chunks = np.array_split(range(len(task_list)), chunks)

            # for i in range(0, len(indices_chunks)):
            #     task_chunks[i]= [task_list[j] for j in indices_chunks[i].tolist()]
            # # print(task_chunks.keys())
            # # print(task_chunks[0])
            # current_chunk = 0

            for i in range(0, len(task_list)):
                current_team.in_queue.put(task_list[i])
            time.sleep(timeout)

            self.logger.info("{} available rules".format(len(sorted_rules)))
            matching_moves = []
            waiting_loops = 0
            # for task in task_chunks[current_chunk]:
            #     current_team.in_queue.put(task)
            # current_chunk = current_chunk + 1

            try:
                # for task in task_chunks[current_chunk]:
                #     current_team.in_queue.put(task)
                # current_chunk = current_chunk + 1
                while True:
                    try:
                        if len(matching_moves) >= self.max_moves:
                            raise FinishedException
                        list_of_moves = current_team.out_queue.get()
                        if list_of_moves != []:
                            #     logging.info("Returned something empty {}".format(list_of_moves))
                            # else:
                            for move_inside_set in list_of_moves:
                                if not move_inside_set.in_list(list_moves = matching_moves, main_layer= main_layer):
                                    matching_moves.append(move_inside_set)
                            # Will need to check here whether it's good or not
                            logging.info("total answer currently is {} for {}".format(matching_moves, self))
                    except Empty:
                        # time.sleep(timeout)
                        # print("nothing new")
                        waiting_loops = waiting_loops + 1

                        # print('Currently empty')
                    # print(len(total_answers))
                    finished_moves = len(matching_moves) >= self.max_moves
                    finished_input_queue = current_team.in_queue.empty()
                    finished_output_queue = current_team.out_queue.empty()
                    if finished_input_queue or finished_output_queue:
                        time.sleep(timeout)
                        print("Extending time because finished input {} or finished output {}".format(finished_input_queue, finished_output_queue))
                    finished_input_queue = current_team.in_queue.empty()
                    # if finished_input_queue:
                    #     print("this should happen at some point")
                    #     for task in task_chunks[current_chunk]:
                    #         print("Currently at chunk {}".format(current_chunk))
                    #         current_team.in_queue.put(task)
                    #     current_chunk = current_chunk + 1
                    finished_output_queue = current_team.out_queue.empty()
                    if finished_moves or (finished_input_queue and finished_output_queue) or (time.time()- start_time > individual_cmp_time):
                        print("stopped because: lengh {}, input {}, output {}, time {}".format(finished_moves, finished_input_queue, finished_output_queue, time.time()- start_time > individual_cmp_time))
                        print("Found {} moves".format(len(matching_moves)))
                        raise FinishedException
                        # break
            except FinishedException:
                self.logger.info("Finished the number of moves")
                self.logger.info("Calculated moves in {} for compound {} - without caching".format(time.time()- start_time, self.InChIKey))
                current_team.flush_team()
                # # print("using killing array")
                #
                # for i in range(len(team)):
                #     killing_array[i] = 1
                # print("other killing streatgy")
                # # in_queue.close()
                # # in_queue.join_thread()
                # # while any([team_player.is_alive() for team_player in team]):
                # #     in_queue.put(POISON_PILL)
                # # in_queue.join()
                # # time.sleep(0.01)
                # for team_player in team:
                #     team_player.join()
                #     print("joining")
                # for team_player in team:
                #     if team_player.is_alive():
                #         self.logger.warning("Some workers are still alive")

        # self.logger.info("For compound {}, matching moves are {}".format(self, matching_moves))
        if use_cache:
            self.logger.info("Caching results for compound {}. The cache contained {}".format(self.InChIKey, len(home_made_cache.keys())))
            home_made_cache[self.InChIKey] = matching_moves
        self.logger.info("Calculated moves in {} for compound {} - without caching".format(time.time()- start_time, self.InChIKey))
        return(matching_moves)

    def apply_transformation_with_move(self, move):
        """
        Used to find applicable rules.
        """
        if move is None:
            self.logger.error("No move is specified")
            return([])
        else:
            return(move.product_list)

    def _get_ECFP(self):
        """
        After comparison of results and speed, decision was made on the following points:
        - Tanimoto is better than Dice
        - Explicit Bit Vectors are better than counted vectors (hashed Morgan fingerprints)
        - Radius 2 Morgan FP will be used
        - Size of bit vector will be 1024 (no difference on similarity on my tests but faster)
        """
        if self.ECFP is None:
            self.ECFP = Chem.AllChem.GetMorganFingerprintAsBitVect(self.rdmol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
        return(self.ECFP)

def __cli():
    """Command line interface for the Compound class"""
    d = "Can check how many rules apply to a compounds and whether a product of interest is present"
    parser = argparse.ArgumentParser(description=d)
    # Logs and saving information
    parser.add_argument("--verbose", help="Default logger is INFO, switch to DEBUG is specified",
                        dest='verbose', action='store_true', default=False)
    parser.add_argument("--log_file", help="Default logger is stderr, switch to log_file if specified",
                        default=None)
    # Compound information
    parser.add_argument("--c_name", help="Compound name. Defaults to None (InchiKey)",
                        default=None)
    # One of the next 2 arguments has to be specified to give a target structure.
    parser.add_argument("--c_smiles", help="Compound smiles", default = None)
    parser.add_argument("--c_inchi", help="Compound inchi", default = None)

    parser.add_argument("--c_smiles_prod", help="Compound smiles looked for in the rules", default = None)
    parser.add_argument("--c_inchi_prod", help="Compound inchi looked for in the rules", default = None)
    # Timeouts on rdkit processes
    parser.add_argument("--fire_timeout", help = "Time allowed for firing one rule on one substrate",
                        type = float, default = 1)
    parser.add_argument("--standardisation_timeout", help = "Time allowed for standardising results from one application",
                        type = float, default = 5)
    parser.add_argument("--expansion_width", help="Maximum number of children", default=5, type=int)
    parser.add_argument("--biological_score_cut_off", default=0.1, type=float)
    parser.add_argument("--substrate_only_score_cut_off", default=0.3, type=float)
    parser.add_argument("--chemical_score_cut_off", default=0.3, type=float)
    parser.add_argument("--diameter", nargs='+',
        help="Diameters to consider", default=[16], type=int)
    parser.add_argument("--EC_filter", nargs='+',
    help="EC numbers to consider for rules", default=None, type=str)
    args = parser.parse_args()
    # Setting up logs
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
        log_writer = open("{}/{}".format(log_file), "w")
        logging.basicConfig(stream=log_writer,
                            level=logging_level,
                            datefmt='%d/%m/%Y %H:%M:%S',
                            format='%(asctime)s -- %(levelname)s -- %(message)s')

    rules, biological_scoring = get_rules_and_score(full_rules_forward_H=full_rules_forward_H,
                                                    full_rules_retro_H=full_rules_retro_H,
                                                    full_rules_forward_no_H=full_rules_forward_no_H,
                                                    full_rules_retro_no_H=full_rules_retro_no_H,
                                                    add_Hs=add_Hs,
                                                    retro=retrosynthesis,
                                                    diameters=args.diameter,
                                                    small=False,
                                                    c_name=args.c_name,
                                                    filtering_EC=args.EC_filter)

    chemical_scoring_configuration = {
        "biological_score_cut_off": args.biological_score_cut_off,
        "substrate_only_score_cut_off": args.substrate_only_score_cut_off,
        "chemical_score_cut_off": args.chemical_score_cut_off}

    biological_score_config = "Using biological cut off at {}".format(args.biological_score_cut_off)
    substrate_only_score_config = "Using substrate only cut off at {}".format(args.substrate_only_score_cut_off)
    chemical_score_config = "Using chemical score cut off at {}".format(args.chemical_score_cut_off)

    root_compound = Compound(csmiles = args.c_smiles,
                            InChI = args.c_inchi,
                            name = args.c_name,
                            max_moves = args.expansion_width,
                            stereo = False,
                            heavy_standardisation = True,
                            fire_timeout = args.fire_timeout,
                            chemical_scoring_configuration = chemical_scoring_configuration,
                            standardisation_timeout = args.standardisation_timeout)
    logging.info("Compound of interest is {}".format(root_compound))

    if args.c_smiles_prod is None and args.c_inchi_prod is None:
        verifying_product = False
    else:
        searched_product = Compound(csmiles = args.c_smiles,
                                InChI = args.c_inchi,
                                name = args.c_name,
                                max_moves = args.expansion_width,
                                stereo = False,
                                heavy_standardisation = True,
                                fire_timeout = args.fire_timeout,
                                chemical_scoring_configuration = chemical_scoring_configuration,
                                standardisation_timeout = args.standardisation_timeout)
        verifying_product = True

    moves = root_compound.obtain_applicable_transformation_with_move(available_rules = rules,
                                            chemical_scorer = "SubandprodChemicalScorer",
                                            biological_scorer = biological_scoring)
    logging.info("Found {} moves with the given settings".format(len(moves)))
    if verifying_product:
        logging.info("Veryfing {} is in the products of {}".format(searched_product, root_compound))
        found = False
        for move in moves:
            if searched_product.in_list(move.product_list):
                logging.info(move.print_all_attributes())
                found = True
        if not found:
            logging.info("Product {} was not found".format(searched_product))

if __name__ == "__main__":
    __cli()
