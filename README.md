
# Monte Carlo Tree Search presentation

The aim of this project is to run a Monte Carlo Tree Search to perform bio-retrosynthesis, compatible with mono-component reaction rules from RetroRules.
The role of each script is detailed below.
Scripts can generally be run from the command line, and have detailed comments for each function.

Detailed docs are in [document_all_options.md](document_all_options.md).

Chemoinformatics choices are detailed in [chemistry_choices.md](chemistry_choices.md).

# Setting conda environment

```bash
conda create --name MCTS python=3.6
source activate MCTS
conda install --channel rdkit rdkit=2019.03.1.0
conda install pytest
conda install pyyaml
```
After cloning this Git repo, please run
```bash
pip install -e .
```
at the root of the package.

Results can be visualised using the stand-alone Scope Viewer available on GitHub at:
```bash
git clone https://github.com/brsynth/scope-viewer.git
```
For using the toxicity calculator:
```bash
conda install scikit-learn=0.19.1
```
For using a database to cache results, you can find it on GitHub:
```bash
conda install pymongo
git clone https://github.com/brsynth/rp3_dcache.git
```
Then run `pip install -e .` at the root of the downloaded package.

# Set-up data files

- You need to specify in the config.py file where you want to store the data (`data_path`).
- Unless specified, it will be stored in the package folder/data
- Organisms will be stored in `[data_path]/organisms` (`data_path` should not end with a /)
while rules in json format with ECFPS (optional) will be in `[data_path]/rules`

Run the following commands:
```bash
python calculate_rule_sets_similarity.py --rule_address_with_H your_rule_address --rule_address_without_H your_rule_address
python calculate_organisms.py
```

# Testing
__Important__: Tests have to be executed in the root folder of the program, which contains the tests folder.
```bash
python change_config.py --use_cache True --add_Hs True
pytest -v
```

# Command line examples

```bash
python change_config.py  \
    --use_cache True
```

```bash
python Tree.py \
    --log_file tree.log \
    --itermax 1000 \
    --expansion_width 10 \
    --time_budget 7200 \
    --max_depth 7 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 0 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubandprodChemicalScorer \
    --virtual_visits 0 \
    --progressive_bias_strategy 0 \
    --diameter 10 12 14 16 \
    --c_name deoxiviolacein \
    --c_inchi "InChI=1S/C20H13N3O2/c24-19-13(18-12-6-2-4-8-16(12)22-20(18)25)9-17(23-19)14-10-21-15-7-3-1-5-11(14)15/h1-10,21H,(H,22,25)(H,23,24)/b18-13+" \
    --folder_to_save deoxi_07_no_H\
    --biological_score_cut_off 0.1 \
    --substrate_only_score_cut_off 0.7 \
    --chemical_score_cut_off 0.7 \
    --minimal_visit_counts 1
```
Expected result from this command is a folder named 'expected_results' containing:
- a text file called 4_results. This allows easy parsing of the fact that 4 results were generated when running the search.
- json files called deoxiviolacein_N.json, with N from 1 to 4: files containing pathways under json format. There is no ordering of those files.
- deoxiviolacein_best.json: contains the best pathway (ie: with the most visits)
- deoxiviolacein_full_scope.json: contains all pathways in the same json file.
- deoxiviolacein_full_tree_for_MCTS.json: contains the Tree under json format. This allows visualisation node by node using our tree viewer available at [https://github.com/brsynth/scope-viewer](https://github.com/brsynth/scope-viewer).
- deoxiviolacein_iteration_N.json: contains the pathway found at iteration N. Folder should contain 4 such files.
- pickles is a folder containing tree_end_search.pkl: the full pickled tree. This can be used for tree extension or analysis.
- results.csv contains the configuration and the results from this configuration (STOP_REASON should be iteration here)
- tree.log contains execution logs.
Except for the json files containing the word 'tree', all json files can be visualised using our scope viewer available at [https://github.com/brsynth/scope-viewer](https://github.com/brsynth/scope-viewer).


# Example for extension and 'normal' search

We expect no result from this search:
```bash
python change_config.py --DB_CACHE True --DB_time 0  --use_cache True
```

```bash
python Tree.py \
    --log_file tree.log \
    --itermax 1000 \
    --expansion_width 10 \
    --time_budget 7200 \
    --max_depth 7 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 0 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubandprodChemicalScorer \
    --virtual_visits 0 --progressive_bias_strategy 0 \
    --diameter 10 12 14 16 \
    --c_name deoxiviolacein \
    --c_inchi "InChI=1S/C20H13N3O2/c24-19-13(18-12-6-2-4-8-16(12)22-20(18)25)9-17(23-19)14-10-21-15-7-3-1-5-11(14)15/h1-10,21H,(H,22,25)(H,23,24)/b18-13+" \
    --folder_to_save test_tree_extension/deoxi_09 \
    --biological_score_cut_off 0.9  \
    --substrate_only_score_cut_off 0.9 \
    --chemical_score_cut_off 0.9 \
    --minimal_visit_counts 1
```

# To rerun from the same Tree with a more tolerant score

The following command will extend the tree by 10 children. What that means is that a node that had 10 children already can have up to 10 other children added. A node that had only 5 can have up to 15 children added (original:10 plus extension:10).
Morevoer, all node scores (visits and values) are reinitialised, as they can change drastically by allowing new rules. Only the structure is conserved, which allows for much faster descent on already expanded nodes.
We expect 1 pathway from this search.
```bash
python change_config.py --DB_CACHE True --DB_time 0  --use_cache True
```

```bash
python Tree.py  \
    --log_file tree.log  \
    --itermax 1000  \
    --expansion_width 10 \
    --time_budget 7200 \
    --max_depth 7 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 0 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubandprodChemicalScorer \
    --virtual_visits 0 \
    --progressive_bias_strategy 0 \
    --diameter 10 12 14 16 \
    --folder_to_save test_tree_extension/deoxi_05 \
    --tree_to_complete end_search \
    --folder_tree_to_complete test_tree_extension/deoxi_09 \
    --biological_score_cut_off 0.1  \
    --substrate_only_score_cut_off 0.5 \
    --chemical_score_cut_off 0.5 \
    --minimal_visit_counts 1
```

# Exploiting the DB

The DB is used as a cache: each time the application of a rule on a compound is run and takes more than DB_time,
it is stored in that database.

- Visualise: [http://localhost:8081][http://localhost:8081]
- Config: in `config.py` file, imported as a module in all scripts

# Supplement finder

The aim of the supplement_finder script is to find potential media supplements that would allow to make other pathways by simple media supplementation. It is currently limited to 1 supplement to avoid combinatorial explosion. It allows for verification of presence in a database of interest (here: Metanetx), previously standardised under the same conditions as the Tree (with or without hydrogens/stereo). 

Please unzip the databases in data/supplement_finder/data before running this script, as well as the search tree in data/supplement_finder/tree_for_testing/TPA/pickles.

Usage:
```bash
python supplement_finder.py --folder_tree_to_complete data/supplement_finder/tree_for_testing/TPA \
--database_address data/supplement_finder/data/metanetx_extracted_inchikeys.json \
--folder_to_save testing_supplement_finder/TPA
```

# Remarks on the config file

- config.py template is located at clean_data/base_config.py and is used to generate a config.py instance used by all jobs.
- config.py is read by all job instances
- editing config.py while jobs are still running will impact not only new jobs but also job that are still running

# Files organisation

- Each file will contain its own class.
- Firing routines are taken from T.D.' reactor package.
- Tests are in a separate folder and have filenames starting with `test` prefix, as well as class names starting with Test.
- Data for the tests is also in the Tests folder.
- DB is taken from T.D.' DBcache.
- Config file contains a number of global parameters needed in various scripts, mostly to decide which features to use (score cut-offs, DB, caching, progressive widening etc)

Important:
- a pickled_data folder needs to be present in the root of the MCTS (for running tests).
- organisms generation (bottom of chemical_compound_state file) should be run once.

# Description of object classes

- Biological scoring: contains obtains that score a biological rule - a move in MCTS (currently, either randomly or from a dictionnary of scores)
- calculate_organisms: run once at set up to extract organisms under correct format.
- calculate_rule_sets_similarity: run once to set up the rules with original substrates and products for similarity calculations. Pickled rules.
- change config: change the config file from command line.
- Chemical Compound state: contains a wrapper around the compound class. It selects the best available moves of a state, checks terminality etc. The basic MCTS object. It also creates the organisms which are states.
- Chemical scoring: utilities for chemical scoring.
- Chemistry choices: documents important chemistry choices made during the project.
- Compound: wrapper around rdkit mol. Basic object containing all transformations, sorting of rules. That is where the chemistry happens. Currently contains archive functions for parallelisation that will be simplified one day.
- compound scoring: at the moment, contains scoring for toxicity bias.
- Config: Allows for configuration of MCTS parameters used throughout the tree running. Also generates the logs for this.
- convert to SBML: converts a JSON file for a pathway to an SBML file.
- data: contains the data necessary to run RetroPath3.0
- expected_results: cotnains expected results from running the first example udner Testing section.
- MCTS_node: has a state as attribute, but also keeps information of father, son, number of visits, rollout results and so on. This is where most algorithmic improvements are encoded.
- Move: contains the move - ie: mostly a wrapper around a chemical rule, with pointing compound and products, as well as scores and RAVE results.
_ organisms: wrapper around organisms used as sinks for the retrosynthetic search.
- pathway: can store pathways, with utils to add compounds and reactions and export as a json. Is also be compatible with scorers.
- pathway_scoring: used to score pathways.
- representation: a class defining how to represent nodes, states and compounds (text vs colors in the Terminal).
- rewarding: a class containing rewards after a rollout.
- Rollout_policies: contains different ways to sample the rollout (random, proportional to the chemical score, the biological score, combination thereof etc)
- rule_set_sexamples: imports rule from csv files. Used for basic examples, does not allows for similarity scoring.
- rule_sets_similarity: imports rules under a similarity compatible format.
- setup: for compatibility with pip.
- supplement_finder: allows finding of supplements in Trees.
- tests: contains tests for installation.
- Tree: performs the tree search: mostly iterating through the MCTS node objects a defined number of times. Utilities for pathways, saving etc. This is the script to use to run the MCTS search.
- tree_viewer: save tree under a viewer compatible format.
- UCT_policies: contains different ways to select the best child (classical, proportional to the chemical score, the biological score, combination thereof etc)
- utilities: contains utilities to apply reaction rules and standardise compounds.

# Rule input formatting

Rules are imported from `rule_sets_similarity` after calculation with `calculate_rule_sets_similarity.py`. The user
can define his own import to replace the default import method.

A rule set as input with `calculate_rule_sets_similarity.py` needs to have the following characteristics...
- be a csv file that is tab delimited.

... and possess the following keys:
- `Rule_ID`
- `Reaction_ID`: ID of the reaction the rule was learned on. Field cannot be empty.
- `Rule_SMARTS`: mono-component reaction SMARTS as described in Duigou et. al., Nucleic Acids Research, 2019.

Optional keys are:

- `Diameter`: diameter around the reaction center. If absent, set to 0.
- `Substrate_ID` and `Substrate_SMILES`: used for chemical scoring. If smiles is given, ID also needs to be provided.
If absent, score will be 1.
- `Product_IDs` and `Product_SMILES`: used for chemical scoring. If smiles is given, ID also needs to be provided.
If absent, score will be 1. Remark: chemical score is disabled if either substrate or product is missing.
- `Rule_SMILES`: A SMILES depiction of the reaction rule. If missing, set to empty string.
- `Score_normalized`: biological score of the reaction. Set to 1 if absent.
- `Reaction_EC_number`: specifies the EC number of the reaction used as template. Set to unspecified if absent.
- `Rule_usage`: possible values are `retro`, `forward` or `both`. Should the reaction be considered for retrosynthesis,
forward usage or both. Set to `both` if unspecified.


# MCTS improvements currently implemented.

- minimal number of visits per child
- using transposition tables (uses too much memory for some reason).
- use RAVE
- progressive bias (biasing exploration by giving initial values to nodes)
- progressive widening: allow number of children roughly proportional to number of visits of the node
- virtual visits: give a number of virtual visits to avoid too much stochasticity in initial evaluations


# Biosensor working example

We expect one result from this search.

```bash
python change_config.py --DB_CACHE True --DB_time 0  --use_cache True --add_Hs True --biosensor True
```

```bash
python Tree.py  \
    --log_file tree.log \
    --itermax 1000  \
    --expansion_width 20 \
    --time_budget 7200 \
    --max_depth 2 \
    --UCT_policy Biochemical_UCT_1 \
    --UCTK 20 \
    --bias_k 0 \
    --k_rave 50 \
    --Rollout_policy Rollout_policy_random_uniform_on_biochemical_multiplication_score \
    --max_rollout 3 \
    --chemical_scoring SubandprodChemicalScorer \
    --virtual_visits 0 \
    --progressive_bias_strategy max_reward  \
    --diameter 10 12 14 16 \
    --c_name pipecolate \
    --c_inchi "InChI=1S/C6H11NO2/c8-6(9)5-3-1-2-4-7-5/h5,7H,1-4H2,(H,8,9)" \
    --folder_to_save pipecolate \
    --EC_filter 1.5.3.7 1.5.3 \
    --biological_score_cut_off 0.1  \
    --substrate_only_score_cut_off 0.9 \
    --chemical_score_cut_off 0.9 \
    --minimal_visit_counts 1
```

# Various remarks

Best move selection:
- using a multiplication of biological score
- and chemical score based on similarity:
    - similarity towards initial compound to order the rules that will be tested
    - similarity also to products for real move ordering after rule has been applied.

Standardisation:
- When loading sink and source compounds, will go through 'heavy' standardisation
- The rest (within the tree search) will go through normal standardisation.
