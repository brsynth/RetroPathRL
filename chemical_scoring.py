"""
Defines the chemical scoring functions
"""

# General utility packages
import random
import itertools  # For all permutations when IDing the best products
import numpy as np  # Allows for simpler calculations on lists
import logging

# Chemistry packages
from rdkit import DataStructs  # For similarity computation

def list_product(combination):
    """
    Calculates the product of all elements from the list.
    Remark: deprecated, use geometric mean instead.
    """
    score = 1
    for tanimoto in combination:
        score = score * tanimoto
    return(score)

def combine_products(product_list, product_list_bis, max_combination = 1000):
    """
    Calculates all possible combinations of products (native and query products).
    Limited to 1000 combinations, knowing that combinations behave as n! with n the number of products.
    """
    combinations = [(x,product_list_bis) for x in itertools.permutations(product_list,len(product_list_bis))]
    if len(combinations) > max_combination:
        combinations = combinations[0:max_combination]
    return(combinations)

def list_geometric_mean(combination):
    """
    Calculates the geometric mean of the array.
    """
    a = np.array(combination)
    return a.prod()**(1.0/len(a))

def tanimoto_product_calc(native_products_ecfp, query_products_ecfp, verbose = False):
    all_scores = []
    if len(native_products_ecfp) != len(query_products_ecfp):
        # Reject rules that do not produce the same number of compounds.
        logging.debug("Rule does not generate the same number of products: native is {} and new is {}".format(len(native_products_ecfp), len(query_products_ecfp)))
        return(-1)
    combinations = combine_products(product_list = native_products_ecfp, product_list_bis = query_products_ecfp)
    score_list = []
    for combination in combinations:
        tanimoto_combination = []
        native, query = combination[0], combination[1]
        for i in range(len(native)):
            tanimoto = DataStructs.cDataStructs.TanimotoSimilarity(native[i], query[i])
            tanimoto_combination.append(tanimoto)
        score_list.append(list_geometric_mean(tanimoto_combination))
    if verbose:
        logging.debug("Score list length is {} and scores {}".format(len(score_list), score_list))
    return(max(score_list))

class ChemicalScoring(object):
    logger = logging.getLogger(__name__)
    """
    Defines Chemical Scorer objects.
    """
    def __init__(self, scoring_function, name = "ChemicalScoring"):
        self.scoring_function = scoring_function
        self.scoring_warning = True
        self.name = name

    def calculate(self, compound, products = None, rule = None, original_substrates_list = None, original_products_list_list = None):
        if original_substrates_list == [None] and (original_products_list_list is None or original_products_list_list == [None]):
            if self.scoring_warning:
                self.scoring_warning = False
                self.logger.warning("Not using chemical scoring for {}. Default is set to 1".format(self.name))
            return(1)
        score, warning = self.scoring_function(compound, products, rule, original_substrates_list, original_products_list_list)
        if not warning is None:
            self.logger.debug(warning)
        return(score)

def pseudo_random(compound, products, rule, original_substrates_list = None, original_products_list_list = None):
    """
    Was used during development.
    """
    warning = None
    if compound.InChIKey == "NBBJYMSMWIIQGU-UHFFFAOYSA-N":
        if rule == "MNXR94682_MNXM821":
            score = 0.99
        elif rule == "MNXR117465_MNXM821":
            score = 0.88
        else:
            score = random.uniform(0,0.75)
    elif compound.InChIKey == "DNIAPMSPPWPWGF-UHFFFAOYSA-N":
        if rule == "MNXR95713_MNXM90191":
            score = 0.80
        elif rule == "MNXR103108_MNXM90191":
            score = 0.76
        else:
            score = random.uniform(0,0.75)
    else:
        score = random.uniform(0,0.75)
    return(score, warning)

def substrate_calculation(compound, products = None, rule = None, original_substrates_list = None, original_products_list_list = None):
    """
    If the original_substrates_list is none, it means chemical scoring is not implemented and scoring should eb neutral: 1 in mutiplication.
    """
    if original_substrates_list is None:
        warning = "Score is set to 1 for cmp {} and rule {}".format(compound, rule)
        tanimoto = 1
    else:
        tanimoto = 0
        for native_substrate in original_substrates_list:
            query_substrate = compound._get_ECFP()
            tanimoto_this = DataStructs.cDataStructs.TanimotoSimilarity(query_substrate, native_substrate)
            warning = None
            tanimoto = max(tanimoto, tanimoto_this)
    return(tanimoto, warning)

def substrate_and_product_calculation(compound, products, rule, original_substrates_list = None, original_products_list_list = None):
    """
    If the original_substrates_list is none, it means chemical scoring is not implemented and scoring should eb neutral: 1 in mutiplication.
    """
    warning = None
    if original_substrates_list is None:
        warning = "Score is set to 1 for cmp {} and rule {}".format(compound, rule)
        tanimoto = 1
        return(tanimoto, warning)
    else:
        tanimoto = 0
        for i in range(len(original_substrates_list)):
            native_substrate = original_substrates_list[i]
            query_substrate = compound._get_ECFP()
            tanimoto_substrate = DataStructs.cDataStructs.TanimotoSimilarity(query_substrate, native_substrate)
            warning = None
            query_products_ecfp = []
            for prod in products:
                query_products_ecfp.append(prod._get_ECFP())
            prod_result = tanimoto_product_calc(original_products_list_list[i], query_products_ecfp, verbose = False)
            if prod_result == -1:
                warning = "Number of product issue with rule {} and products {}".format(rule, products)
            tanimoto = max(tanimoto, tanimoto_substrate * prod_result)
        return(tanimoto, warning)

def constant_scorer(compound, products, rule, original_substrates_list = None, original_products_list_list = None):
    warning = None
    return(1, warning)

RandomChemicalScorer = ChemicalScoring(scoring_function = pseudo_random, name = "RandomChemicalScorer")
SubstrateChemicalScorer = ChemicalScoring(scoring_function = substrate_calculation, name = "SubstrateChemicalScorer")
SubandprodChemicalScorer = ChemicalScoring(scoring_function = substrate_and_product_calculation, name = "SubandprodChemicalScorer")
ConstantChemicalScorer = ChemicalScoring(scoring_function = constant_scorer, name = "ConstantChemicalScorer")
# Chemical scoring utilities. Taken from similarity.
