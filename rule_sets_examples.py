"""
Contains the rules examples that will be used throughout the tests.
The aim is to
"""

import logging
import csv
import os

rule_10_subset_address = "{}/tests/data/rules_r10_subset.tsv".format(os.path.dirname(__file__))
applicable_rules_10_dict = {}
with open(rule_10_subset_address, "r") as csv_file:
    fieldnames = ["Rule_ID", "Reaction_ID", "Diameter", "Direction", "Rule_order", "Rule_SMARTS", "Substrate_ID", "Substrate_SMILES", "Product_IDs", "Product_SMILES", "Rule_SMILES", "Rule_SMARTS_lite"]
    csv_reader = csv.DictReader(csv_file, delimiter = '\t', fieldnames = fieldnames)
    next(csv_reader)  # skip first line
    for element in csv_reader:
        applicable_rules_10_dict[element["Rule_ID"]] = {"Rule_SMARTS": element["Rule_SMARTS"],
                                                        "biological_score": 1,
                                                        "EC_number": ["EC: None"],
                                                        "Rule_SMILES": element["Rule_SMILES"]}


rule_2_subset_address = "{}/tests/data/rules_r2_subset.tsv".format(os.path.dirname(__file__))
applicable_rules_2_dict = {}
with open(rule_2_subset_address, "r") as csv_file:
    fieldnames = ["Rule_ID", "Reaction_ID", "Diameter", "Direction", "Rule_order", "Rule_SMARTS", "Substrate_ID", "Substrate_SMILES", "Product_IDs", "Product_SMILES", "Rule_SMILES", "Rule_SMARTS_lite"]
    csv_reader = csv.DictReader(csv_file, delimiter = '\t', fieldnames = fieldnames)
    next(csv_reader)  # skip first line
    for element in csv_reader:
        applicable_rules_2_dict[element["Rule_ID"]] =  {"Rule_SMARTS": element["Rule_SMARTS"],
                                                        "biological_score": 1,
                                                        "EC_number": ["EC: None"],
                                                        "Rule_SMILES": element["Rule_SMILES"]}


rule_mixed_subset_address = "{}/tests/data/rules_mixed_subset.tsv".format(os.path.dirname(__file__))
applicable_rules_mixed_dict = {}
with open(rule_mixed_subset_address, "r") as csv_file:
    fieldnames = ["Rule_ID", "Reaction_ID", "Diameter", "Direction", "Rule_order", "Rule_SMARTS", "Substrate_ID", "Substrate_SMILES", "Product_IDs", "Product_SMILES", "Rule_SMILES", "Rule_SMARTS_lite"]
    csv_reader = csv.DictReader(csv_file, delimiter = '\t', fieldnames = fieldnames)
    next(csv_reader)  # skip first line
    for element in csv_reader:
        applicable_rules_mixed_dict[element["Rule_ID"]] = {"Rule_SMARTS": element["Rule_SMARTS"],
                                                        "biological_score": 1,
                                                        "EC_number": ["EC: None"],
                                                        "Rule_SMILES": element["Rule_SMILES"]}
