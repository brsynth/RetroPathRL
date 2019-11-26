# Documentation

The aim of this file is to document all options available to run the MCTS and where to find them.
More details are in the attached paper at https://doi.org/10.1101/800474 , especially in the Appendix.

### Global configuration options

- where: in the config.py file
- how: either by modying the config.py file by hand or by running the change_config.py with its argparser (recommanded)

- DB_CACHE: uses the MongoDB cache when activated
- DB_REPLACE: replaces data in the Mongo DB cache when activated
- DB_time: time cut_off for loading in the DB: stored only if above the cut-off, otherwise the rule is applied by Python
- use_cache: dictionnary for caching results within the script. Highly recommanded.
- retrosynthesis: performs a retrosynthetic search; biosensor: performs a biosensor search
Both cannot be activated at the same time.
Main difference is how the state is considered sucessful: all compounds have to be found for retrosynthesis, and only 1 for biosensors
- add_Hs: explicit hydrogens. Recommanded at false for faster calculations.
- use_transpositions and transposition depth: not stable. Allow for sharing of information between nodes with the same chemical state but at different places in the tree, as done in doi:10.1007/BF03192151.

### Tree search configuration:

- stop_at_first_result: stops once a signle pathway is found.
- c_name, c_smile,s c_inchi: information on the chemical compound of interest
- fire_timeout, standardisation_timeout: time allowed for firing a rule/standardising a compound
- organism_name: which model to use for production of compounds
- complementary_sink: csv file containing compounds to add to the sink. If organism is None, is the full sink.
- representation: how to print results in logs
- itermax: maximum number of iterations allowed for running the Tree search
- parallel: not possible to use at the moment due to workaround for RDKit rule application. Aimed at parallelising rollouts.
- expansion_width: maximum number of children per node
- time budget: time allowed for running the tree search. The search will stop at the end of the iteration that exceeds this allotted time
- max_depth: maximum depth of the Tree (also the maximum number of pathway steps)
- minimal_visit_counts: Minimal number of times a node has to be rolled out before his brothers can be expanded
- UCT_policy: define the UCT policy to use, ie the way to rank the children of a node. Allows various bias, and scoring considerations.
- UCTK: the constant defining the exploration/exploitation parameter in the UCT formula
- bias_k: if progressive bias is used, define the weight of the progressive bias in the UCT formula
- k_rave: if RAVE is used, how to weight the RAVE. Roughly for visits below this value RAVE values lead the UCT and above, rollout values lead.
- use_RAVE: moves have scores each time they are used throughout the Tree, adapting RAVE (Rapid Action Value Estimation) principle to the whole tree and not just rollouts.
- penalty: penalty when no compound of the state belongs to the organism 
- full_state_reward: reward when all compounds of the state belongs to the organism 
- pathway_scoring: how to score a pathway when it is found.
- Rollout_policy: how to select moves for the Rollout: randomly, wieghting by which scores. Many options.
- max_rollout: maximum length of the rollout (it also stops when max_depth is reached)
- chemical_scoring: chose the way to chemically score reactions (considering only substrates or both substrates and products). Possibility to use ConstantChemicalScorer which always returns 1.
- biological_score_cut_off: cuts off with biological score at the specified level 
- substrate_only_score_cut_off: cuts off with substrate similarity only score BEFORE applying the rule at the specified level 
- chemical_score_cut_off: cuts off with specified chemical score AFTER applying the rule at the specified level 
- virtual_visits: start nodes at virtual_visits values, to avoid stochasticity at initial simulations; used to avoid having to much variability at initial Monte Carlo simulations.
- progressive_bias_strategy: policy for the progressive bias (untested)
- progressive widening: add a child to nodes visited more than len(nodes)^2 (untested)
- diameter: Speficy the diameters (as list) to use
- EC_filter: allow only certain EC subclasses
- small: development archive
- seed: for reproducibility
- tree_to_complete: if restarting the search from another tree.



