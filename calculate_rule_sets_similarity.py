"""
Calculates pickled rules for rule dataset, including compounds ECFP for similarity.
The aim is to load rules taking into account:
- original substrates and products with their ECFPs
- normalised biological score
- pickle the rules so that this is calculated only once at the first use of RP3.
"""

# General utilities
import logging
import csv
import pickle
import os
import sys
import argparse

# Chemistry packages
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
# Specific RP3 packages
from utilities.reactor.Core import RuleBurnerCore
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError

from config import *


class IncorrectFormatting(Exception):
    """Class for incorrect loading of rules."""

    def __init__(self, message):
        self.message = message


def __run__(rule_address_with_H = None, rule_address_without_H = None, rm_stereo = True):


    all_ecfps = {}

    rule_10_subset_address = "{}/tests/data/rules_r10_subset.tsv".format(os.path.dirname(os.path.abspath(__file__)))
    applicable_rules_10_dict_sim = {}
    with open(rule_10_subset_address, "r") as csv_file:
        fieldnames = ["Rule_ID", "Reaction_ID", "Diameter", "Direction", "Rule_order", "Rule_SMARTS", "Substrate_ID", "Substrate_SMILES", "Product_IDs", "Product_SMILES", "Rule_SMILES", "Rule_SMARTS_lite"]
        csv_reader = csv.DictReader(csv_file, delimiter = '\t', fieldnames = fieldnames)
        next(csv_reader)  # skip first line
        for element in csv_reader:
            applicable_rules_10_dict_sim[element["Rule_ID"]] = {"Rule_SMARTS": element["Rule_SMARTS"],
                                                            "biological_score": 1,
                                                            "EC_number": ["EC: None"],
                                                            "Rule_SMILES": element["Rule_SMILES"]}
            # Obtaining substrate ECFP:
            if element["Substrate_ID"] in all_ecfps.keys():
                ECFP = all_ecfps[element["Substrate_ID"]]
            else:
                rd_mol = Chem.MolFromSmiles(element["Substrate_SMILES"], sanitize=False)  # Important: Sanitize = False
                rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                all_ecfps[element["Substrate_ID"]] = ECFP
            applicable_rules_10_dict_sim[element["Rule_ID"]]["substrate_ECFP"] = ECFP
            product_names = element["Product_IDs"].split(".")
            product_smiles = element["Product_SMILES"].split(".")
            applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"] = []
            for i in range(len(product_smiles)):
                if product_names[i] in all_ecfps.keys():
                    ECFP = all_ecfps[product_names[i]]
                else:
                    rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                    rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                    ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                    all_ecfps[product_names[i]] = ECFP
                applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"].append(ECFP)


    rule_mixed_subset_address = "{}/tests/data/rules_mixed_subset.tsv".format(os.path.dirname(os.path.abspath(__file__)))
    applicable_rules_mixed_dict_sim = {}
    with open(rule_mixed_subset_address, "r") as csv_file:
        fieldnames = ["Rule_ID", "Reaction_ID", "Diameter", "Direction", "Rule_order", "Rule_SMARTS", "Substrate_ID", "Substrate_SMILES", "Product_IDs", "Product_SMILES", "Rule_SMILES", "Rule_SMARTS_lite"]
        csv_reader = csv.DictReader(csv_file, delimiter = '\t', fieldnames = fieldnames)
        next(csv_reader)  # skip first line
        for element in csv_reader:
            applicable_rules_mixed_dict_sim[element["Rule_ID"]] = {"Rule_SMARTS": element["Rule_SMARTS"],
                                                            "biological_score": 1,
                                                            "EC_number": ["EC: None"],
                                                            "Rule_SMILES": element["Rule_SMILES"]}
            if element["Substrate_ID"] in all_ecfps.keys():
                ECFP = all_ecfps[element["Substrate_ID"]]
            else:
                rd_mol = Chem.MolFromSmiles(element["Substrate_SMILES"], sanitize=False)  # Important: Sanitize = False
                rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                all_ecfps[element["Substrate_ID"]] = ECFP
            applicable_rules_mixed_dict_sim[element["Rule_ID"]]["substrate_ECFP"] = ECFP
            product_names = element["Product_IDs"].split(".")
            product_smiles = element["Product_SMILES"].split(".")
            applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"] = []
            for i in range(len(product_smiles)):
                if product_names[i] in all_ecfps.keys():
                    ECFP = all_ecfps[product_names[i]]
                else:
                    rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                    rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                    ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                    all_ecfps[product_names[i]] = ECFP
                applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"].append(ECFP)

    def merge_rule_characteristics(current_characteristics, new_characteristics, rule_id):
        """
        The aim of this function is to merge new_characteristics of a rule into current_characteristics.
        It will be called when parsing real datasets.
        This is because rule IDs are based on hashed SMARTS, meaning rules from differetn original reaction can have the smae ID.
        It checks the
        - used diameters
        - susbtrate and products ECFPs
        - Reaction_ID
        - Biological score
        """
        assert current_characteristics["Rule_SMARTS"] == new_characteristics["Rule_SMARTS"]
        assert current_characteristics["Diameter"] == new_characteristics["Diameter"]
        # If compounds and products are already in list, do not add them
        common_substrate = new_characteristics["Substrate_ID"][0] in current_characteristics["Substrate_ID"]
        common_products = new_characteristics["Product_IDs"][0] in current_characteristics["Product_IDs"]
        if common_substrate and common_products:
            logging.info("Hashed rule {} versions ({} and {}) are based on the same substrates and products".format(rule_id,
                                                                                new_characteristics["Reaction_ID"][0],
                                                                                current_characteristics["Reaction_ID"]))
        else:
            current_characteristics["Substrate_ID"].append(new_characteristics["Substrate_ID"][0])
            current_characteristics["substrate_ECFP"].append(new_characteristics["substrate_ECFP"][0])
            current_characteristics["Product_IDs"].append(new_characteristics["Product_IDs"][0])
            current_characteristics["products_ECFP"].append(new_characteristics["products_ECFP"][0])
        # Biological score:
        if current_characteristics["biological_score"] != new_characteristics["biological_score"]:
            logging.info("Different biological scores for {}".format(rule_id))
            current_characteristics["biological_score"] = max(current_characteristics["biological_score"], new_characteristics["biological_score"])
        # Merging EC numbers
        for EC in new_characteristics["EC_number"]:
            if EC not in current_characteristics["EC_number"]:
                current_characteristics["EC_number"].append(EC)
        # Merging Original Identifiers
        if new_characteristics["Reaction_ID"][0] not in current_characteristics["Reaction_ID"]:
            current_characteristics["Reaction_ID"].append(new_characteristics["Reaction_ID"][0])
        return(current_characteristics)



    logging.info("-------------------Currently processing rules with Hs--------------------------")
    all_ecfps = {}
    full_rules_retro_H = {}
    full_rules_forward_H = {}

    if rule_address_with_H is None:
        logging.warning("No rule address with H was given. Moving on to the next dataset")
    else:
        with open(rule_address_with_H, "r") as csv_file:
            csv_reader = csv.DictReader(csv_file, delimiter = '\t')
            next(csv_reader)
            # Following booleans are for logging missing columns only once.
            # They are useful when using rule formats without as much information as RetroRules.
            first_log_biological_None = True
            first_log_EC_None = True
            first_log_rule_usage_None = True
            first_log_diameter_None = True
            first_log_rule_smiles_None= True
            first_log_chemistry_sub_None= True
            first_log_chemistry_sub_name_None = True
            first_log_chemistry_prod_None= True
            first_log_chemistry_prod_name_None = True
            for element in csv_reader:
                try:
                    rule_id = element["Rule_ID"]
                except KeyError:
                    try:
                        rule_id = element["# Rule_ID"]
                    except KeyError:
                        message = "Rule_ID is mandatory. If you do not have such a column, please create it"
                        raise IncorrectFormatting(message)
                try:
                    reac_id = element["Reaction_ID"]
                except KeyError:
                    message = "Reaction_ID is mandatory. If you do not have such a column, please create it"
                    raise IncorrectFormatting(message)
                try:
                    try:
                        biological_score = float(element["Score_normalized"])
                    except ValueError:
                        logging.warning("No biological score for {}. Is set to 1".format(rule_id))
                except KeyError:
                    if first_log_biological_None:
                        first_log_biological_None = False
                        logging.warning("No column Score_normalized. Default biological score is set to 1")
                    biological_score = 1

                # Calculate ECFPs - substrate
                try:
                    sub_smiles = element["Substrate_SMILES"]
                    try:
                        sub_id = element["Substrate_ID"]
                        # ID and smiles present, calculate chemistry.
                        if sub_id in all_ecfps.keys():
                            sub_ECFP = all_ecfps[sub_id]
                        else:
                            rd_mol = Chem.MolFromSmiles(sub_smiles, sanitize=False)  # Important: Sanitize = False
                            rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                            sub_ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                            all_ecfps[sub_id] = sub_ECFP
                    except KeyError:
                        message = "Substrate_ID is mandatory is you specify Substrate_SMILES \n"
                        message_bis = "If you do not have such a column, please create it. You can use InChIKey if you do not have a proper ID"
                        raise IncorrectFormatting(message + message_bis)
                except KeyError:
                    if first_log_chemistry_sub_None:
                        first_log_chemistry_sub_None = False
                        logging.warning("No column Substrate_SMILES. Set to None \n Chemical scoring will not be available")
                    sub_smiles = None
                    sub_ECFP = None
                    try:
                        sub_id = element["Substrate_ID"]
                    except KeyError:
                        if first_log_chemistry_sub_name_None:
                            first_log_chemistry_sub_name_None = False
                            logging.warning("No column Substrate_ID. Set to Unspecified")
                        sub_id = "Unspecified"
                # # Calculate ECFPs - products
                try:
                    product_smiles = element["Product_SMILES"].split(".")
                    try:
                        product_ids = element["Product_IDs"].split(".")
                        assert len(product_smiles) == len(product_ids)
                        # All information is avalaible for ECFP calculation for products
                        products_ECFPs = []
                        for i in range(len(product_smiles)):
                            if product_ids[i] in all_ecfps.keys():
                                prod_ECFP = all_ecfps[product_ids[i]]
                            else:
                                rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                                rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=rm_stereo)
                                prod_ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                                all_ecfps[product_ids[i]] = prod_ECFP
                            products_ECFPs.append(prod_ECFP)
                        product_ids = [set(element['Product_IDs'].split('.'))]
                    except KeyError:
                        message = "Product_IDs is mandatory is you specify Product_SMILES \n"
                        message_bis = "If you do not have such a column, please create it. You can use InChIKeys if you do not have proper IDs"
                        raise IncorrectFormatting(message + message_bis)
                    except AssertionError:
                        message = "Product IDs and structures do not have the same len for {} ({})".format(product_names, rule_id)
                        raise IncorrectFormatting(message)
                except KeyError:
                    if first_log_chemistry_prod_None:
                        first_log_chemistry_prod_None = False
                        logging.warning("No column product_smiles. Set to None \n Chemical scoring will not be available")
                    product_smiles = None
                    products_ECFPs = None
                    try:
                        product_ids = element["Product_IDs"]
                    except KeyError:
                        if first_log_chemistry_prod_name_None:
                            first_log_chemistry_prod_name_None = False
                            logging.warning("No column Product_IDs. Set to Unspecified")
                        product_ids = "Unspecified"
                # EC number:
                try:
                    EC_number = element["Reaction_EC_number"].split(",")
                except KeyError:
                    if first_log_EC_None:
                        first_log_EC_None = False
                        logging.warning("No column Reaction_EC_number. Default EC is set to Unspecified")
                    EC_number = "Unspecified"

                # Rule characteristics
                # Keep smart, keep original ids as lists
                try:
                    diameter = element["Diameter"]
                except KeyError:
                    if first_log_diameter_None:
                        first_log_diameter_None = False
                        logging.warning("No column Diameter. Default diameter is set to 0")
                    diameter = 0
                # Load rule smarts
                try:
                    Rule_SMARTS = element["Rule_SMARTS"]
                except KeyError:
                    message = "Rule_SMARTS is mandatory. If you do not have such a column, please create it"
                    raise IncorrectFormatting(message)
                try:
                    Rule_SMILES = element["Rule_SMILES"]
                except KeyError:
                    if first_log_rule_smiles_None:
                        first_log_rule_smiles_None = False
                        logging.warning("No column Rule_SMILES. Default Rule_SMILES is set to empty string")
                    Rule_SMILES = ""
                rule_characteristics_extended = {"Rule_SMARTS": Rule_SMARTS,
                                                "biological_score": biological_score,
                                                "EC_number": EC_number,
                                                "substrate_ECFP": [sub_ECFP],
                                                "products_ECFP": [products_ECFPs],
                                                "Reaction_ID" : [reac_id],
                                                "Substrate_ID" : [sub_id],
                                                "Product_IDs" : product_ids,
                                                "Diameter" : diameter,
                                                "Rule_SMILES" : Rule_SMILES
                                                }
                try:
                    usage = element["Rule_usage"]
                except KeyError:
                    if first_log_rule_usage_None:
                        logging.warning("No column Rule_usage. Default usage is set to both directions")
                        first_log_rule_usage_None = False
                    usage = "both"

                if usage == "forward":
                    if rule_id in full_rules_forward_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_forward_H[rule_id],
                                                                    new_characteristics = rule_characteristics_extended,
                                                                    rule_id = rule_id)
                        full_rules_forward_H[rule_id] = new_characteristics
                    else:
                        full_rules_forward_H[rule_id] = rule_characteristics_extended
                elif usage == "retro":
                    if rule_id in full_rules_retro_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_retro_H[rule_id],
                                                                new_characteristics = rule_characteristics_extended,
                                                                rule_id = rule_id)
                        full_rules_retro_H[rule_id] = new_characteristics
                    else:
                        full_rules_retro_H[rule_id] = rule_characteristics_extended
                elif usage == "both":
                    # Adding to retro
                    if rule_id in full_rules_retro_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_retro_H[rule_id],
                                                                new_characteristics = rule_characteristics_extended,
                                                                rule_id = rule_id)
                        full_rules_retro_H[rule_id] = new_characteristics
                    else:
                        full_rules_retro_H[rule_id] = rule_characteristics_extended
                    # Adding to forward
                    if rule_id in full_rules_forward_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_forward_H[rule_id],
                                                                    new_characteristics = rule_characteristics_extended,
                                                                    rule_id = rule_id)
                        full_rules_forward_H[rule_id] = new_characteristics
                    else:
                        full_rules_forward_H[rule_id] = rule_characteristics_extended
                else:
                    raise NotImplementedError

        useful_characteristics = ["Rule_SMARTS",
                                  "biological_score",
                                  "EC_number",
                                  "substrate_ECFP",
                                  "products_ECFP",
                                  "Reaction_ID",
                                  "Diameter",
                                  "Rule_SMILES"]

        for rid, rule_char in full_rules_forward_H.items():
            full_rules_forward_H[rid] = dict((attribute, rule_char[attribute]) for attribute in useful_characteristics)

        for rid, rule_char in full_rules_retro_H.items():
            full_rules_retro_H[rid] = dict((attribute, rule_char[attribute]) for attribute in useful_characteristics)

        with open("{}/{}.pkl".format(folder_to_save, "full_rules_forward_H"), "wb") as pickle_handler:
            pickle.dump(full_rules_forward_H, pickle_handler)

        with open("{}/{}.pkl".format(folder_to_save, "full_rules_retro_H"), "wb") as pickle_handler:
            pickle.dump(full_rules_retro_H, pickle_handler)


    logging.info("-------------------Currently processing rules without Hs--------------------------")

    all_ecfps_no_H = {}
    full_rules_retro_no_H = {}
    full_rules_forward_no_H = {}

    if rule_address_without_H is None:
        logging.warning("No rule address without H was given.")
    else:
        with open(rule_address_without_H, "r") as csv_file:
            csv_reader = csv.DictReader(csv_file, delimiter = '\t')
            next(csv_reader)
            first_log_biological_None = True
            first_log_EC_None = True
            first_log_rule_usage_None = True
            first_log_diameter_None = True
            first_log_rule_smiles_None = True
            first_log_chemistry_sub_None= True
            first_log_chemistry_sub_name_None = True
            first_log_chemistry_prod_None= True
            first_log_chemistry_prod_name_None = True

            for element in csv_reader:
                try:
                    rule_id = element["Rule_ID"]
                except KeyError:
                    try:
                        rule_id = element["# Rule_ID"]
                    except KeyError:
                        message = "Rule_ID is mandatory. If you do not have such a column, please create it"
                        raise IncorrectFormatting(message)
                try:
                    try:
                        biological_score = float(element["Score_normalized"])
                    except ValueError:
                        logging.warning("No biological score for {}. Is set to 1".format(rule_id))
                        biological_score = 1
                except KeyError:
                    if first_log_biological_None:
                        first_log_biological_None = False
                        logging.warning("No column Score_normalized. Default biological score is set to 1")
                    biological_score = 1

                # Calculate ECFPs - substrate
                try:
                    sub_smiles = element["Substrate_SMILES"]
                    try:
                        sub_id = element["Substrate_ID"]
                        # ID and smiles present, calculate chemistry.
                        if sub_id in all_ecfps_no_H.keys():
                            sub_ECFP = all_ecfps_no_H[sub_id]
                        else:
                            rd_mol = Chem.MolFromSmiles(sub_smiles, sanitize=False)  # Important: Sanitize = False
                            rd_mol = standardize_chemical(rd_mol, add_hs=False, rm_stereo=rm_stereo)
                            sub_ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                            all_ecfps_no_H[sub_id] = sub_ECFP
                    except KeyError:
                        message = "Substrate_ID is mandatory is you specify Substrate_SMILES \n"
                        message_bis = "If you do not have such a column, please create it. You can use InChIKey if you do not have a proper ID"
                        raise IncorrectFormatting(message + message_bis)
                except KeyError:
                    if first_log_chemistry_sub_None:
                        first_log_chemistry_sub_None = False
                        logging.warning("No column Substrate_SMILES. Set to None. \n Chemical scoring will not be available")
                    sub_smiles = None
                    sub_ECFP = None
                    try:
                        sub_id = element["Substrate_ID"]
                    except KeyError:
                        if first_log_chemistry_sub_name_None:
                            first_log_chemistry_sub_name_None = False
                            logging.warning("No column Substrate_ID. Set to Unspecified")
                        sub_id = "Unspecified"
                ## Calculate ECFP:
                try:
                    product_smiles = element["Product_SMILES"].split(".")
                    try:
                        product_ids = element["Product_IDs"].split(".")
                        assert len(product_smiles) == len(product_ids)
                        # All information is avalaible for ECFP calculation for products
                        products_ECFPs = []
                        for i in range(len(product_smiles)):
                            if product_ids[i] in all_ecfps_no_H.keys():
                                prod_ECFP = all_ecfps_no_H[product_ids[i]]
                            else:
                                rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                                rd_mol = standardize_chemical(rd_mol, add_hs=False, rm_stereo=rm_stereo)
                                prod_ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = not rm_stereo)
                                all_ecfps[product_ids[i]] = prod_ECFP
                            products_ECFPs.append(prod_ECFP)
                        product_ids = [set(element['Product_IDs'].split('.'))]
                    except KeyError:
                        message = "Product_IDs is mandatory is you specify Product_SMILES \n"
                        message_bis = "If you do not have such a column, please create it. You can use InChIKeys if you do not have proper IDs"
                        raise IncorrectFormatting(message + message_bis)
                    except AssertionError:
                        message = "Product IDs and structures do not have the same len for {} ({})".format(product_names, rule_id)
                        raise IncorrectFormatting(message)
                except KeyError:
                    if first_log_chemistry_prod_None:
                        first_log_chemistry_prod_None = False
                        logging.warning("No column product_smiles. Set to None \n Chemical scoring will not be available")
                    product_smiles = None
                    products_ECFPs = None
                    try:
                        product_ids = element["Product_IDs"]
                    except KeyError:
                        if first_log_chemistry_prod_name_None:
                            first_log_chemistry_prod_name_None = False
                            logging.warning("No column Product_IDs. Set to Unspecified")
                        product_ids = "Unspecified"
                # EC number:
                try:
                    EC_number = element["Reaction_EC_number"].split(",")
                except KeyError:
                    if first_log_EC_None:
                        first_log_EC_None = False
                        logging.warning("No column Ec_number. Default EC is set to Unspecified")
                    EC_number = "Unspecified"

                # Rule characteristics
                try:
                    diameter = element["Diameter"]
                except KeyError:
                    if first_log_diameter_None:
                        first_log_diameter_None = False
                        logging.warning("No column Diameter. Default diameter is set to 0")
                    diameter = 0
                try:
                    reac_id = element["Reaction_ID"]
                except KeyError:
                    message = "Reaction_ID is mandatory. If you do not have such a column, please create it"
                    raise IncorrectFormatting(message)
                try:
                    Rule_SMARTS = element["Rule_SMARTS"]
                except KeyError:
                    message = "Rule_SMARTS is mandatory. If you do not have such a column, please create it"
                    raise IncorrectFormatting(message)
                try:
                    Rule_SMILES = element["Rule_SMILES"]
                except KeyError:
                    if first_log_rule_smiles_None:
                        first_log_rule_smiles_None = False
                        logging.warning("No column Rule_SMILES. Default Rule_SMILES is set to empty string")
                    Rule_SMILES = ""


                rule_characteristics_extended = {"Rule_SMARTS": Rule_SMARTS,
                                                "biological_score": biological_score,
                                                "EC_number": EC_number,
                                                "substrate_ECFP": [sub_ECFP],
                                                "products_ECFP": [products_ECFPs],
                                                "Reaction_ID" : [reac_id],
                                                "Substrate_ID" : [sub_id],
                                                "Product_IDs" : product_ids,
                                                "Diameter" : diameter,
                                                "Rule_SMILES" : Rule_SMILES
                                                }
                try:
                    usage = element["Rule_usage"]
                except KeyError:
                    if first_log_rule_usage_None:
                        logging.warning("No column Rule_usage. Default usage is set to both directions")
                        first_log_rule_usage_None = False
                    usage = "both"

                if usage == "forward":
                    if rule_id in full_rules_forward_no_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_forward_no_H[rule_id],
                                                                    new_characteristics = rule_characteristics_extended,
                                                                    rule_id = rule_id)
                        full_rules_forward_no_H[rule_id] = new_characteristics
                    else:
                        full_rules_forward_no_H[rule_id] = rule_characteristics_extended
                elif usage == "retro":
                    if rule_id in full_rules_retro_no_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_retro_no_H[rule_id],
                                                                new_characteristics = rule_characteristics_extended,
                                                                rule_id = rule_id)
                        full_rules_retro_no_H[rule_id] = new_characteristics
                    else:
                        full_rules_retro_no_H[rule_id] = rule_characteristics_extended
                elif usage == "both":
                    # Adding to retro
                    if rule_id in full_rules_retro_no_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_retro_no_H[rule_id],
                                                                new_characteristics = rule_characteristics_extended,
                                                                rule_id = rule_id)
                        full_rules_retro_no_H[rule_id] = new_characteristics
                    else:
                        full_rules_retro_no_H[rule_id] = rule_characteristics_extended
                    # Adding to forward
                    if rule_id in full_rules_forward_no_H.keys():
                        new_characteristics = merge_rule_characteristics(current_characteristics = full_rules_forward_no_H[rule_id],
                                                                    new_characteristics = rule_characteristics_extended,
                                                                    rule_id = rule_id)
                        full_rules_forward_no_H[rule_id] = new_characteristics
                    else:
                        full_rules_forward_no_H[rule_id] = rule_characteristics_extended

                else:
                    raise NotImplementedError

        useful_characteristics = ["Rule_SMARTS",
                                  "biological_score",
                                  "EC_number",
                                  "substrate_ECFP",
                                  "products_ECFP",
                                  "Reaction_ID",
                                  "Diameter",
                                  "Rule_SMILES"]

        for rid, rule_char in full_rules_forward_no_H.items():
            full_rules_forward_no_H[rid] = dict((attribute, rule_char[attribute]) for attribute in useful_characteristics)

        for rid, rule_char in full_rules_retro_no_H.items():
            full_rules_retro_no_H[rid] = dict((attribute, rule_char[attribute]) for attribute in useful_characteristics)

        with open("{}/{}.pkl".format(folder_to_save, "full_rules_forward_no_H"), "wb") as pickle_handler:
            pickle.dump(full_rules_forward_no_H, pickle_handler)

        with open("{}/{}.pkl".format(folder_to_save, "full_rules_retro_no_H"), "wb") as pickle_handler:
            pickle.dump(full_rules_retro_no_H, pickle_handler)


if __name__ == "__main__":
    d = "Formatting rules in a RP3 compatible format"
    parser = argparse.ArgumentParser(description=d)
    # Logs and saving information
    parser.add_argument("--terminal", help="Default logger is logs_rules_set_up, switch to terminal if specified",
                        action='store_true', default=False)
    parser.add_argument("--rule_address_with_H", default = None,
                        help="Rule file with hydrogen that will be used.")
    parser.add_argument("--rule_address_without_H", default = None,
                        help="Rule file without hydrogen that will be used")
    parser.add_argument("--stereo", help="By default, removing stereo. This parameter allows for stereo enabling",
                        action='store_true', default=False)

    args = parser.parse_args()

    folder_to_save = "{}/{}".format(data_path, "rules")
    if not os.path.exists(folder_to_save):
        os.mkdir(folder_to_save)

    if args.terminal is True:
        logging.basicConfig(
                stream = sys.stderr,
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
    else:
        logging.basicConfig(
                stream = open("{}/{}.log".format(folder_to_save, "logs_rules_set_up"), "w"),
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
        print("By default, logs are saved in {}/logs_rules_set_up.log. Please use --terminal to redirect to sys.stderr".format(folder_to_save))
    try:
        if args.stereo:
            logging.warning("Stereo is enabled")
        else:
            logging.info("Removing stereo")
        __run__(rule_address_with_H = args.rule_address_with_H, rule_address_without_H = args.rule_address_without_H, rm_stereo = not args.stereo)
    except IncorrectFormatting as e:
        logging.error(e.message)
        print("Inccorect termination. See logs for more details")
