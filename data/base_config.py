"""
the aim of this file is to store configuration parameters, notably for the DB.
It replaces what I previously wanted to define as 'global'
"""
try:
    from rp3_dcache.Manager import Manager  # In house module
    from rp3_dcache.Utils import make_document_id, as_document, rdmols_from_document
    dcache_installed = True
except ModuleNotFoundError:
    dcache_installed = False
import logging
import os

# Files and addresses configurations - should not be modified:
global data_path
data_path = "{}/data".format(os.path.dirname(os.path.abspath(__file__)))

global add_Hs
add_Hs = True
hydrogen_config = "Using explicit hydrogens : {}".format(add_Hs)

# Database for storing results configuration
global DB_CACHE
global DB_REPLACE
DB_CACHE = False and dcache_installed
DB_REPLACE = False and dcache_installed
DB_time = 0
if DB_CACHE:
    global CACHE_MGR
    if add_Hs:
        CACHE_MGR = Manager(replace=DB_REPLACE, collection = "results_with_H")
    else:
        CACHE_MGR = Manager(replace=DB_REPLACE, collection = "results_without_H")
    CACHE_MGR.connect()
    DB_config = "Setting the DB from config file: Installed package: {}. Using cache DB: {}; Replacing results: {}".format(dcache_installed, DB_CACHE, DB_REPLACE)
elif dcache_installed:
    DB_config = "Setting the DB from config file: Installed package: {}. Using cache DB: {}; Replacing results: {}".format(dcache_installed, DB_CACHE, DB_REPLACE)
else:
    DB_config = "Setting the DB from config file: Installed package: {}".format(dcache_installed)

# Mode for using RP3: retrosynthesis or biosensor. QSAR might be implemented one day.
global retrosynthesis
global biosensor
retrosynthesis = True
biosensor = False
tree_mode_config = "Using retrosynthesis: {} - using biosensor {}".format(retrosynthesis, biosensor)

# Configuring local cache. Could be replaced by a proper caching system one day.
global home_made_cache
home_made_cache = {}

global use_cache
use_cache = True

cache_config = "Initialising an empty cache: {}; Using it: {}".format(home_made_cache, use_cache)

# MCTS parameters for configuration

global transposition_table
global use_transpositions
global use_transpositions_depth

transposition_table = {}
use_transpositions = False
use_transpositions_depth = False

transposition_table_config = "Using transposition tables: {}. With depth: {}".format(use_transpositions, use_transpositions_depth)

# For toxicity, using log(IC50) as penaly when below 0.
global use_toxicity
try:
    import sklearn
    from sklearn.neural_network import MLPRegressor
    sklearn_here = True
except ModuleNotFoundError:
    toxicity_config = "Toxicity will not be enabled because sklearn is not installed"
    sklearn_here = False
use_toxicity = False
use_toxicity = use_toxicity and sklearn_here
