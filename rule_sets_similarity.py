"""
Contains the rules examples that will be used throughout the tests.
The aim is to
"""

import logging
import csv
import pickle
import os
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
# To standardise compounds before calculating the said ECFPs
from utilities.reactor.Core import RuleBurnerCore
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError
from config import *

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
            rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=True)
            ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
            all_ecfps[element["Substrate_ID"]] = ECFP
        applicable_rules_10_dict_sim[element["Rule_ID"]]["substrate_ECFP"] = [ECFP]
        product_names = element["Product_IDs"].split(".")
        product_smiles = element["Product_SMILES"].split(".")
        applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"] = []
        for i in range(len(product_smiles)):
            if product_names[i] in all_ecfps.keys():
                ECFP = all_ecfps[product_names[i]]
            else:
                rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=True)
                ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
                all_ecfps[product_names[i]] = ECFP
            applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"].append(ECFP)
        applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"] = [applicable_rules_10_dict_sim[element["Rule_ID"]]["products_ECFP"]]


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
            rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=True)
            ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
            all_ecfps[element["Substrate_ID"]] = ECFP
        applicable_rules_mixed_dict_sim[element["Rule_ID"]]["substrate_ECFP"] = [ECFP]
        product_names = element["Product_IDs"].split(".")
        product_smiles = element["Product_SMILES"].split(".")
        applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"] = []
        for i in range(len(product_smiles)):
            if product_names[i] in all_ecfps.keys():
                ECFP = all_ecfps[product_names[i]]
            else:
                rd_mol = Chem.MolFromSmiles(product_smiles[i], sanitize=False)  # Important: Sanitize = False
                rd_mol = standardize_chemical(rd_mol, add_hs=True, rm_stereo=True)
                ECFP = AllChem.GetMorganFingerprintAsBitVect(rd_mol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
                all_ecfps[product_names[i]] = ECFP
            applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"].append(ECFP)
        applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"] = [applicable_rules_mixed_dict_sim[element["Rule_ID"]]["products_ECFP"]]


folder_for_data = "{}/{}".format(data_path, "rules")

with open("{}/{}.pkl".format(folder_for_data, "full_rules_forward_H"), "rb") as pickle_handler:
    full_rules_forward_H = pickle.load(pickle_handler)

with open("{}/{}.pkl".format(folder_for_data, "full_rules_retro_H"), "rb") as pickle_handler:
    full_rules_retro_H = pickle.load(pickle_handler)

with open("{}/{}.pkl".format(folder_for_data, "full_rules_forward_no_H"), "rb") as pickle_handler:
    full_rules_forward_no_H = pickle.load(pickle_handler)

with open("{}/{}.pkl".format(folder_for_data, "full_rules_retro_no_H"), "rb") as pickle_handler:
    full_rules_retro_no_H = pickle.load(pickle_handler)

# Information below is used when small is activated.
TPA_keys = [
'RR-01-ff15c6f282fb-10-F_10', 'RR-01-7467598facd5-10-F_10', 'RR-01-eced43368870-10-F_10',
'RR-01-4be04021e193-10-F_10', 'RR-01-855e527e8dac-10-F_10', 'RR-01-be6c9619b09d-10-F_10',
'RR-01-6d06e27aef2e-10-F_10', 'RR-01-ba2eb2715516-10-F_10', 'RR-01-3791e719f852-10-F_10',
'RR-01-c7ebbe91808e-10-F_10', 'RR-01-2ebdcd7010bd-10-F_10', 'RR-01-f765514b772f-10-F_10',
'RR-01-3545d6ed7403-10-F_10', 'RR-01-2c7b5bc4bd50-10-F_10', 'RR-01-56602d66630b-10-F_10',
'RR-01-0b9242207244-10-F_10', 'RR-01-6ae5e8c91fd8-10-F_10', 'RR-01-73e66a22ce11-10-F_10',
'RR-01-a397b384391a-10-F_10', 'RR-01-bdfadc58b5a9-10-F_10', 'RR-01-398ea892cfda-10-F_10',
'RR-01-d55dd578e73e-10-F_10', 'RR-01-4c79f297821e-10-F_10', 'RR-01-ec985dcd998f-10-F_10',
'RR-01-f88421e31b7f-10-F_10', 'RR-01-c05e8a6eaa8e-10-F_10', 'RR-01-439c65246142-10-F_10',
'RR-01-915d8ffeafc7-10-F_10', 'RR-01-1312d40b0eca-10-F_10', 'RR-01-db16fa6490be-10-F_10',
'RR-01-211fa91ad9f6-10-F_10', 'RR-01-d3a021cc85d1-10-F_10', 'RR-01-f348f5dd6f53-10-F_10',
'RR-01-6402887eee8c-12-F_12', 'RR-01-1b0e63da9155-10-F_10', 'RR-01-eef619ff8017-10-F_10']

TPA_16_keys = [
    'RR-02-0b92422072445abd-16-F', 'RR-02-a397b384391a4069-16-F', 'RR-02-ff15c6f282fb8c05-10-F', 'RR-02-7467598facd5dc87-10-F', 'RR-02-eced43368870d3d5-10-F', 'RR-02-4be04021e1933793-10-F', 'RR-02-855e527e8dacdf66-10-F', 'RR-02-be6c9619b09df557-16-F', 'RR-02-6d06e27aef2ecfca-10-F', 'RR-02-ba2eb2715516342a-10-F', 'RR-02-2ebdcd7010bd7d71-10-F', 'RR-02-3545d6ed74034d56-16-F', 'RR-02-2c7b5bc4bd503c71-16-F', 'RR-02-56602d66630bcfc3-16-F', 'RR-02-6ae5e8c91fd8d8af-16-F', 'RR-02-73e66a22ce113ef0-16-F', 'RR-02-bdfadc58b5a9d059-16-F', 'RR-02-398ea892cfda59a1-16-F', 'RR-02-d55dd578e73e0d77-16-F', 'RR-02-4c79f297821e39c8-16-F', 'RR-02-ec985dcd998f1ad5-16-F', 'RR-02-f88421e31b7f840d-10-F', 'RR-02-c05e8a6eaa8e5d7b-10-F', 'RR-02-439c65246142bb8d-16-F', 'RR-02-eef619ff80171bb9-10-F', 'RR-02-1b0e63da91559964-10-F', 'RR-02-915d8ffeafc759e3-16-F', 'RR-02-1312d40b0eca0500-10-F', 'RR-02-db16fa6490be4edc-16-F', 'RR-02-211fa91ad9f6118e-16-F', 'RR-02-d3a021cc85d1df38-10-F', 'RR-02-f348f5dd6f534894-10-F', 'RR-02-6402887eee8ccf7f-12-F'
]

TPA_keys_archive = ['MNXR10030_MNXM1095', 'MNXR11026_MNXM217', 'MNXR14800_MNXM1095', 'MNXR14815_MNXM217', 'MNXR14815_MNXM26', 'MNXR15002_MNXM1403', 'MNXR15002_MNXM438', 'MNXR17596_MNXM217', 'MNXR17596_MNXM60', 'MNXR17596_MNXM2734', 'MNXR7145_MNXM438', 'MNXR7162_MNXM29', 'MNXR7162_MNXM371', 'MNXR7162_MNXM4645', 'MNXR8092_MNXM1095', 'MNXR8092_MNXM987', 'MNXR8092_MNXM3766', 'MNXR8092_MNXM438', 'MNXR81819_MNXM1132', 'MNXR81819_MNXM240', 'MNXR81819_MNXM92', 'MNXR84236_MNXM217', 'MNXR84236_MNXM92', 'MNXR84236_MNXM240', 'MNXR87313_MNXM1095', 'MNXR87313_MNXM210', 'MNXR93560_MNXM210', 'MNXR93560_MNXM4645', 'MNXR9766_MNXM217', 'MNXR4531_MNXM3500', 'MNXR242_MNXM4077', 'MNXR16010_MNXM3331']

violacein_keys = [
'RR-01-73b9e2dc9ddf-12-F_12', 'RR-01-40f5dd9f4609-12-F_12', 'RR-01-3ab86ef33d83-12-F_12',
'RR-01-3d2f047d8c5d-12-F_12', 'RR-01-6cd1048f02ac-12-F_12', 'RR-01-79b5cb776128-12-F_12',
'RR-01-ddd8913a0977-12-F_12', 'RR-01-c8f4eb4b319f-12-F_12', 'RR-01-f46c2cc731e3-12-F_12',
'RR-01-1dacab2272d6-12-F_12', 'RR-01-96b3271dca32-12-F_12', 'RR-01-dba62271635d-12-F_12',
'RR-01-6f464e721366-12-F_12', 'RR-01-22668a073a24-12-F_12', 'RR-01-7bf8c911187c-12-F_12',
'RR-01-eb021cb9b45f-12-F_12', 'RR-01-ea0941b7e233-12-F_12', 'RR-01-b8a979b1e0af-12-F_12',
'RR-01-b416c220e769-12-F_12', 'RR-01-b761dec4bc0b-12-F_12']

violacein_16_keys = [
    'RR-02-73b9e2dc9ddf8783-12-F', 'RR-02-40f5dd9f46097fca-12-F', 'RR-02-ddd8913a09771dc1-16-F',
    'RR-02-c8f4eb4b319f1271-16-F', 'RR-02-f46c2cc731e3f3e6-12-F', 'RR-02-1dacab2272d62122-12-F',
    'RR-02-96b3271dca322f14-12-F', 'RR-02-dba62271635df6d8-12-F', 'RR-02-6f464e721366959e-16-F',
    'RR-02-22668a073a248f56-16-F', 'RR-02-7bf8c911187cfdb8-16-F', 'RR-02-eb021cb9b45f1c77-16-F',
    'RR-02-ea0941b7e2330bb7-12-F', 'RR-02-b8a979b1e0aff9db-12-F', 'RR-02-b416c220e7694280-16-F',
    'RR-02-b761dec4bc0b6d69-16-F', 'RR-02-efb2ae8a1e477ce7-02-F'
]

violacein_keys_archive = ['MNXR2613_MNXM1378', 'MNXR62938_MNXM97338', 'MNXR62939_MNXM97621', 'MNXR62940_MNXM88916', 'MNXR62941_MNXM97503', 'MNXR62942_MNXM98206', 'MNXR62943_MNXM98911', 'MNXR62944_MNXM98919', 'MNXR76570_MNXM10640', 'MNXR85789_MNXM316', 'MNXR85789_MNXM61', 'MNXR85789_MNXM7216', 'MNXR87814_MNXM10640', 'MNXR87814_MNXM97268']

styrene_keys = [
'RR-01-3545d6ed7403-16-F_16', 'RR-01-2c7b5bc4bd50-16-F_16', 'RR-01-90b72f2f7271-16-F_16',
'RR-01-86e6fb5af7e7-16-F_16', 'RR-01-d55dd578e73e-16-F_16', 'RR-01-4c79f297821e-16-F_16',
'RR-01-9ea719dd7a38-16-F_16', 'RR-01-f71a10986373-16-F_16', 'RR-01-c63f6198a7e6-16-F_16',
'RR-01-90207d9c857e-16-F_16', 'RR-01-5f2baf1968dc-16-F_16', 'RR-01-051398bdbd09-16-F_16',
'RR-01-285481018aad-16-F_16', 'RR-01-ec985dcd998f-16-F_16', 'RR-01-a5663bfcd012-16-F_16',
'RR-01-e90f89cdba81-16-F_16', 'RR-01-439c65246142-16-F_16', 'RR-01-438622e469b9-16-F_16',
'RR-01-0ce86297b8a0-16-F_16', 'RR-01-7f682d0b3544-16-F_16', 'RR-01-caf68bc1dcd9-16-F_16',
'RR-01-80cb6058ffaf-16-F_16', 'RR-01-fd0c6fc209da-16-F_16', 'RR-01-4f0ee4666374-16-F_16',
'RR-01-0da874351cfe-16-F_16', 'RR-01-2ab0796934b2-16-F_16', 'RR-01-c3ec74eeb9cb-16-F_16',
'RR-01-80c1be875426-16-F_16', 'RR-01-0d46c69dc939-16-F_16', 'RR-01-564c599b8a7d-16-F_16',
'RR-01-0ba3992e57de-16-F_16', 'RR-01-e7e3ea6fcb2e-16-F_16']

styrene_16_keys = [
    'RR-02-3545d6ed74034d56-16-F', 'RR-02-2c7b5bc4bd503c71-16-F', 'RR-02-d55dd578e73e0d77-16-F', 'RR-02-4c79f297821e39c8-16-F', 'RR-02-051398bdbd090ac0-16-F', 'RR-02-285481018aad5970-16-F', 'RR-02-ec985dcd998f1ad5-16-F', 'RR-02-439c65246142bb8d-16-F', 'RR-02-c3ec74eeb9cbf299-16-F', 'RR-02-80c1be875426e6d6-16-F', 'RR-02-c63f6198a7e6f821-16-F', 'RR-02-438622e469b97ba2-16-F', 'RR-02-caf68bc1dcd9ce0f-16-F', 'RR-02-80cb6058ffaff5b7-16-F', 'RR-02-0ba3992e57dea8f9-16-F', 'RR-02-e7e3ea6fcb2e454e-16-F', 'RR-02-90207d9c857ede9f-16-F', 'RR-02-4f0ee46663746d19-16-F', 'RR-02-0d46c69dc9394d1e-16-F', 'RR-02-564c599b8a7d6009-16-F', 'RR-02-90b72f2f72710436-16-F', 'RR-02-86e6fb5af7e720f6-16-F', 'RR-02-9ea719dd7a385419-16-F', 'RR-02-f71a109863731e17-16-F', 'RR-02-5f2baf1968dcb2a8-16-F', 'RR-02-a5663bfcd012ad53-16-F', 'RR-02-e90f89cdba81cdae-16-F', 'RR-02-0ce86297b8a08bd6-16-F', 'RR-02-7f682d0b3544ff05-16-F', 'RR-02-fd0c6fc209da098a-16-F', 'RR-02-0da874351cfe02c9-16-F', 'RR-02-2ab0796934b264ab-16-F'
]

styrene_keys_archive = ['MNXR10030_MNXM1095', 'MNXR14800_MNXM1095', 'MNXR14852_MNXM4534', 'MNXR15002_MNXM438', 'MNXR20707_MNXM8128', 'MNXR23023_MNXM2434', 'MNXR23044_MNXM438', 'MNXR6957_MNXM2108', 'MNXR7145_MNXM438', 'MNXR7765_MNXM2787', 'MNXR8092_MNXM438', 'MNXR8581_MNXM2034', 'MNXR86322_MNXM890', 'MNXR87207_MNXM987', 'MNXR87313_MNXM1095', 'MNXR93681_MNXM438', 'MNXR94080_MNXM2434']

pinocembrin_keys = [
'RR-01-ccf71c47bb9d-12-F_12', 'RR-01-0c8072624504-12-F_12', 'RR-01-1c5839a0ab53-12-F_12',
'RR-01-f8a09ffeaff2-12-F_12', 'RR-01-89f31ace6816-12-F_12', 'RR-01-3545d6ed7403-12-F_12',
'RR-01-2c7b5bc4bd50-12-F_12', 'RR-01-84f2d75c7006-12-F_12', 'RR-01-d8260cc37483-12-F_12',
'RR-01-5d0e22c36664-12-F_12', 'RR-01-619f5488435d-12-F_12', 'RR-01-ace320d40808-12-F_12',
'RR-01-1eef64c0807c-12-F_12', 'RR-01-90c8db21b14a-12-F_12', 'RR-01-7a232ac9d78d-12-F_12',
'RR-01-fcfba06ea22e-12-F_12', 'RR-01-0ad3ebab6842-12-F_12', 'RR-01-e00ad30f92ef-12-F_12',
'RR-01-0bdebd001527-12-F_12', 'RR-01-be7eaf23b012-12-F_12', 'RR-01-713c9f21bb7a-12-F_12',
'RR-01-fa5526899164-12-F_12', 'RR-01-698e74ec6c01-12-F_12', 'RR-01-455f8c3cb357-12-F_12',
'RR-01-2feefffd698e-12-F_12', 'RR-01-080bc4abd404-12-F_12', 'RR-01-8cc405a7228a-12-F_12',
'RR-01-94c48631eb63-12-F_12', 'RR-01-4a2d24769cd2-12-F_12', 'RR-01-fe4dd69d1d04-12-F_12',
'RR-01-d9836ab8784e-12-F_12', 'RR-01-8cee5d752624-12-F_12', 'RR-01-d6480cb0e6c1-12-F_12',
'RR-01-9b87c1cf20bf-12-F_12', 'RR-01-814cbdaa7721-12-F_12', 'RR-01-e65449ee4a2b-12-F_12',
'RR-01-57b8cbef05c8-12-F_12', 'RR-01-21e16356a12d-12-F_12', 'RR-01-a36521afcc76-12-F_12']

pinocembrin_16_keys = [
    'RR-02-3545d6ed74034d56-16-F', 'RR-02-2c7b5bc4bd503c71-16-F', 'RR-02-90c8db21b14a0cdd-16-F', 'RR-02-7a232ac9d78ddca1-16-F', 'RR-02-8cc405a7228a62f2-16-F', 'RR-02-94c48631eb63e101-16-F', 'RR-02-9b87c1cf20bf1721-16-F', 'RR-02-ccf71c47bb9df145-12-F', 'RR-02-0c8072624504c7b6-12-F', 'RR-02-1c5839a0ab53e203-12-F', 'RR-02-f8a09ffeaff2866a-12-F', 'RR-02-89f31ace6816319f-12-F', 'RR-02-84f2d75c7006398e-16-F', 'RR-02-d8260cc3748381a8-16-F', 'RR-02-5d0e22c36664b1a9-12-F', 'RR-02-619f5488435db48a-12-F', 'RR-02-ace320d40808fb1e-12-F', 'RR-02-1eef64c0807c929d-12-F', 'RR-02-fcfba06ea22ecd80-16-F', 'RR-02-0ad3ebab6842f767-16-F', 'RR-02-e00ad30f92ef189d-12-F', 'RR-02-0bdebd00152790b2-12-F', 'RR-02-be7eaf23b012eaa1-12-F', 'RR-02-fa55268991642d87-12-F', 'RR-02-698e74ec6c012ee8-12-F', 'RR-02-2feefffd698ec847-12-F', 'RR-02-4a2d24769cd20003-16-F', 'RR-02-fe4dd69d1d04fe31-16-F', 'RR-02-d9836ab8784e7c31-16-F', 'RR-02-8cee5d752624b8ec-16-F', 'RR-02-d6480cb0e6c1ddfa-12-F', 'RR-02-814cbdaa77217477-12-F', 'RR-02-e65449ee4a2bbe2e-12-F', 'RR-02-57b8cbef05c8057c-12-F', 'RR-02-21e16356a12d9e8e-16-F', 'RR-02-a36521afcc766305-16-F'
]
pinocembrin_keys_archive = ['MNXR7145_MNXM438', 'MNXR93681_MNXM438', 'MNXR93681_MNXM97', 'MNXR1041_MNXM443', 'MNXR2251_MNXM406', 'MNXR227_MNXM264', 'MNXR14993_MNXM953', 'MNXR60189_MNXM3188', 'MNXR84871_MNXM1208', 'MNXR85701_MNXM1794', 'MNXR85702_MNXM5784', 'MNXR60602_MNXM4497', 'MNXR60602_MNXM45001', 'MNXR73989_MNXM34231', 'MNXR73989_MNXM6269', 'MNXR84948_MNXM556', 'MNXR85459_MNXM1025', 'MNXR85703_MNXM2889']

naringenin_keys = [
'RR-01-ccf71c47bb9d-12-F_12', 'RR-01-0c8072624504-12-F_12', 'RR-01-1c5839a0ab53-12-F_12',
'RR-01-f8a09ffeaff2-12-F_12', 'RR-01-89f31ace6816-12-F_12', 'RR-01-3545d6ed7403-12-F_12',
'RR-01-2c7b5bc4bd50-12-F_12', 'RR-01-84f2d75c7006-12-F_12', 'RR-01-d8260cc37483-12-F_12',
'RR-01-5d0e22c36664-12-F_12', 'RR-01-619f5488435d-12-F_12', 'RR-01-2fd0f61d4287-12-F_12',
'RR-01-be4344437152-12-F_12', 'RR-01-ace320d40808-12-F_12', 'RR-01-1eef64c0807c-12-F_12',
'RR-01-90c8db21b14a-12-F_12', 'RR-01-7a232ac9d78d-12-F_12', 'RR-01-fcfba06ea22e-12-F_12',
'RR-01-0ad3ebab6842-12-F_12', 'RR-01-e00ad30f92ef-12-F_12', 'RR-01-0bdebd001527-12-F_12',
'RR-01-be7eaf23b012-12-F_12', 'RR-01-713c9f21bb7a-12-F_12', 'RR-01-fa5526899164-12-F_12',
'RR-01-698e74ec6c01-12-F_12', 'RR-01-455f8c3cb357-12-F_12', 'RR-01-2feefffd698e-12-F_12',
'RR-01-080bc4abd404-12-F_12', 'RR-01-8cc405a7228a-12-F_12', 'RR-01-94c48631eb63-12-F_12',
'RR-01-4a2d24769cd2-12-F_12', 'RR-01-fe4dd69d1d04-12-F_12', 'RR-01-a5d6105dddf3-12-F_12',
'RR-01-3dc4bc5586d2-12-F_12', 'RR-01-d9836ab8784e-12-F_12', 'RR-01-8cee5d752624-12-F_12',
'RR-01-19bd510ccb82-12-F_12', 'RR-01-8c06c2f87518-12-F_12', 'RR-01-d6480cb0e6c1-12-F_12',
'RR-01-9b87c1cf20bf-12-F_12', 'RR-01-814cbdaa7721-12-F_12', 'RR-01-e65449ee4a2b-12-F_12',
'RR-01-57b8cbef05c8-12-F_12', 'RR-01-21e16356a12d-12-F_12', 'RR-01-a36521afcc76-12-F_12',
'RR-01-27b915ea8462-49-F_-2', 'RR-01-d7542c98f32f-49-F_-2', 'RR-01-bc668c943bcb-49-F_-2',
'RR-01-ebab689ac2f6-49-F_-2']

naringenin_16_keys = [
'RR-02-3545d6ed74034d56-16-F', 'RR-02-2c7b5bc4bd503c71-16-F', 'RR-02-90c8db21b14a0cdd-16-F', 'RR-02-7a232ac9d78ddca1-16-F', 'RR-02-8cc405a7228a62f2-16-F', 'RR-02-94c48631eb63e101-16-F', 'RR-02-9b87c1cf20bf1721-16-F', 'RR-02-ccf71c47bb9df145-12-F', 'RR-02-0c8072624504c7b6-12-F', 'RR-02-1c5839a0ab53e203-12-F', 'RR-02-f8a09ffeaff2866a-12-F', 'RR-02-89f31ace6816319f-12-F', 'RR-02-84f2d75c7006398e-16-F', 'RR-02-d8260cc3748381a8-16-F', 'RR-02-5d0e22c36664b1a9-12-F', 'RR-02-619f5488435db48a-12-F', 'RR-02-2fd0f61d4287301d-12-F', 'RR-02-be4344437152172f-12-F', 'RR-02-ace320d40808fb1e-12-F', 'RR-02-1eef64c0807c929d-12-F', 'RR-02-d7542c98f32fdee5-16-F', 'RR-02-fcfba06ea22ecd80-16-F', 'RR-02-0ad3ebab6842f767-16-F', 'RR-02-e00ad30f92ef189d-12-F', 'RR-02-0bdebd00152790b2-12-F', 'RR-02-be7eaf23b012eaa1-12-F', 'RR-02-fa55268991642d87-12-F', 'RR-02-698e74ec6c012ee8-12-F', 'RR-02-2feefffd698ec847-12-F', 'RR-02-4a2d24769cd20003-16-F', 'RR-02-fe4dd69d1d04fe31-16-F', 'RR-02-a5d6105dddf324f4-12-F', 'RR-02-3dc4bc5586d209c1-12-F', 'RR-02-d9836ab8784e7c31-16-F', 'RR-02-8cee5d752624b8ec-16-F', 'RR-02-bc668c943bcb0817-16-F', 'RR-02-19bd510ccb8276fb-12-F', 'RR-02-8c06c2f87518a9ea-12-F', 'RR-02-d6480cb0e6c1ddfa-12-F', 'RR-02-814cbdaa77217477-12-F', 'RR-02-e65449ee4a2bbe2e-12-F', 'RR-02-57b8cbef05c8057c-12-F', 'RR-02-ebab689ac2f69cfa-16-F', 'RR-02-27b915ea8462b7d9-16-F', 'RR-02-21e16356a12d9e8e-16-F', 'RR-02-a36521afcc766305-16-F'
]

naringenin_keys_archive = ['MNXR15003_MNXM987', 'MNXR17733_MNXM5901', 'MNXR18415_MNXM264', 'MNXR227_MNXM264', 'MNXR7111_MNXM505', 'MNXR7145_MNXM438', 'MNXR84871_MNXM1208', 'MNXR84948_MNXM556', 'MNXR85459_MNXM1025', 'MNXR8874_MNXM505', 'MNXR93681_MNXM438', 'MNXR93685_MNXM505']

ethylene_keys = [
'RR-01-b618db3f535a-16-F_16', 'RR-01-5bed8bcb898a-16-F_16', 'RR-01-c89a1f53ed77-16-F_16',
'RR-01-8226f8c194da-16-F_16', 'RR-01-4d6f36c61098-16-F_16', 'RR-01-a4e51784e182-16-F_16',
'RR-01-1ff6dce7b962-16-F_16', 'RR-01-f11ba70e708f-16-F_16', 'RR-01-340f1860db05-16-F_16',
'RR-01-6f335b33a0a1-16-F_16', 'RR-01-466420eb1f9f-16-F_16', 'RR-01-6059f029037e-16-F_16',
'RR-01-3c44623b991f-16-F_16', 'RR-01-76e374dd28d8-16-F_16', 'RR-01-3f7be56e3f9f-16-F_16',
'RR-01-b052b51fefa9-16-F_16', 'RR-01-8c1d80d26269-16-F_16', 'RR-01-ce0036d18be1-16-F_16',
'RR-01-d386eb98d569-16-F_16', 'RR-01-456711777040-16-F_16', 'RR-01-6122c11077a9-16-F_16',
'RR-01-38bbc744ca2a-16-F_16', 'RR-01-0933bc9cdbc0-16-F_16', 'RR-01-f75e0e113579-16-F_16',
'RR-01-d71acaf8e1fd-16-F_16', 'RR-01-47ecc93f52bc-16-F_16', 'RR-01-234e5a367f8d-16-F_16',
'RR-01-19b55fc4074e-16-F_16', 'RR-01-993c9253f78c-16-F_16', 'RR-01-53f7dd114286-16-F_16',
'RR-01-4d30d6315e8e-16-F_16', 'RR-01-b8b0271c3111-16-F_16', 'RR-01-71834150ed73-16-F_16',
'RR-01-6454a14b930b-16-F_16', 'RR-01-b87dd759d81b-16-F_16', 'RR-01-cae9afc75fed-16-F_16',
'RR-01-b9ad3b92215e-16-F_16', 'RR-01-f2fb1c5efcb7-16-F_16', 'RR-01-a15388df9112-16-F_16',
'RR-01-afeda5dd9e00-16-F_16', 'RR-01-ba1333507ae6-16-F_16', 'RR-01-c8661b3fd89f-16-F_16',
'RR-01-f2a89d1fc110-16-F_16', 'RR-01-61a59a875941-16-F_16', 'RR-01-69bd577ca2f4-16-F_16',
'RR-01-bcd763d5f2e0-16-F_16', 'RR-01-f59041f27577-16-F_16', 'RR-01-860446453e33-16-F_16',
'RR-01-10c14f73d7c8-16-F_16', 'RR-01-b76dd8e8d3e1-16-F_16', 'RR-01-910859fe14f1-16-F_16',
'RR-01-b8ac6a83aeab-16-F_16', 'RR-01-c8d7128110cc-16-F_16', 'RR-01-3f6b8212b11c-16-F_16',
'RR-01-17f237b35b9c-16-F_16', 'RR-01-ad1fcf135fe7-16-F_16', 'RR-01-056cd46e7fbe-16-F_16',
'RR-01-347b66837a93-16-F_16', 'RR-01-2fa9f69671f3-16-F_16', 'RR-01-1c5e2dd6014d-16-F_16',
'RR-01-fe69aaf88d7c-16-F_16', 'RR-01-b8b5ef05dd6f-16-F_16', 'RR-01-b4b4e7eff745-16-F_16',
'RR-01-ba05f56e7531-16-F_16', 'RR-01-c469571bf5dd-16-F_16', 'RR-01-6be59bc48287-16-F_16',
'RR-01-7508e4dad853-16-F_16', 'RR-01-0a1283e064f6-16-F_16', 'RR-01-0590ce05c110-16-F_16',
'RR-01-06663cb0780a-16-F_16', 'RR-01-746e3ff5aa9c-16-F_16', 'RR-01-a609d1ef085d-16-F_16',
'RR-01-91911ecdb81c-16-F_16', 'RR-01-094f718f9fbf-16-F_16', 'RR-01-65c0adfb547c-16-F_16',
'RR-01-c3ec74eeb9cb-16-F_16', 'RR-01-80c1be875426-16-F_16', 'RR-01-8818dea65ff3-16-F_16',
'RR-01-dc579a8e1f8c-16-F_16', 'RR-01-af8d6139571f-16-F_16', 'RR-01-e810bb4aa2b8-16-F_16',
'RR-01-4aa0893bf4c8-16-F_16', 'RR-01-ee7ccc9f6296-16-F_16', 'RR-01-fec67f6f4132-16-F_16',
'RR-01-a3dab3e47518-16-F_16', 'RR-01-b387337dc94e-16-F_16', 'RR-01-1c65c477fa18-16-F_16',
'RR-01-532930687155-16-F_16', 'RR-01-60c1aab00243-16-F_16', 'RR-01-052f180261ba-16-F_16',
'RR-01-2355de6b0a7f-16-F_16', 'RR-01-73f4f0b40f22-16-F_16', 'RR-01-757af66bae93-16-F_16',
'RR-01-ab535f7af6f9-16-F_16', 'RR-01-d3ff173f19d9-16-F_16', 'RR-01-92ea06368d8a-16-F_16',
'RR-01-5d459efbf36a-16-F_16', 'RR-01-0b06353b386c-16-F_16', 'RR-01-c1baf484736a-16-F_16',
'RR-01-eddc4c82bc7f-16-F_16', 'RR-01-42464e412b6a-16-F_16', 'RR-01-d56f2788cec9-16-F_16',
'RR-01-a27c10121672-16-F_16', 'RR-01-14c9ef7ac8a2-16-F_16', 'RR-01-aa26c97ac100-16-F_16',
'RR-01-5549dac72492-16-F_16', 'RR-01-4004b7ea811f-16-F_16', 'RR-01-2704be5eea79-16-F_16',
'RR-01-7a86b7668b9e-16-F_16', 'RR-01-7e4de265f24b-16-F_16', 'RR-01-2ab116f84d68-16-F_16',
'RR-01-e14ce618aee8-16-F_16', 'RR-01-469888aee72f-16-F_16', 'RR-01-c2c66cedf094-16-F_16',
'RR-01-0be9b0efc9e1-16-F_16', 'RR-01-3d3a456f054a-16-F_16', 'RR-01-ac6d35e02a43-16-F_16',
'RR-01-3b31aefba79f-16-F_16', 'RR-01-1e8e39bd8644-16-F_16', 'RR-01-085db4297707-16-F_16',
'RR-01-c5375175a751-16-F_16', 'RR-01-721183f2c545-16-F_16', 'RR-01-031bd961221f-16-F_16',
'RR-01-40baa9f9769c-16-F_16', 'RR-01-dd377c6c9f41-16-F_16', 'RR-01-d60e6efb94fb-16-F_16',
'RR-01-40f28277bfc5-16-F_16', 'RR-01-03b5810baa31-16-F_16', 'RR-01-cc9f5226009b-16-F_16',
'RR-01-289a74d34301-16-F_16', 'RR-01-fbcc8cfd08a6-16-F_16', 'RR-01-2068030457f1-16-F_16',
'RR-01-8b4324fcfd2b-16-F_16', 'RR-01-b62f37c6a404-16-F_16', 'RR-01-4fe900904d0a-16-F_16',
'RR-01-3ba0719fa57a-16-F_16', 'RR-01-7c7ab7ac6b98-16-F_16', 'RR-01-6fec6966f1f2-16-F_16',
'RR-01-36bb1f194ab4-16-F_16', 'RR-01-643b041ce09c-16-F_16', 'RR-01-db0cdb3811be-16-F_16',
'RR-01-e68f663510b3-16-F_16', 'RR-01-bfbf4361736f-16-F_16']

ethylene_16_keys = [
    'RR-02-ba1333507ae6006e-16-F', 'RR-02-c8661b3fd89f1ba4-16-F', 'RR-02-056cd46e7fbe267e-16-F', 'RR-02-347b66837a93582a-16-F', 'RR-02-0a1283e064f61a88-16-F', 'RR-02-0590ce05c110cf1a-16-F', 'RR-02-f11ba70e708f7f3e-16-F', 'RR-02-6059f029037e4eb4-16-F', 'RR-02-3c44623b991f1a2b-16-F', 'RR-02-0933bc9cdbc054bc-16-F', 'RR-02-06663cb0780a094c-16-F', 'RR-02-746e3ff5aa9c9fe9-16-F', 'RR-02-4aa0893bf4c817a5-16-F', 'RR-02-ee7ccc9f62967de0-16-F', 'RR-02-8226f8c194dae8a7-16-F', 'RR-02-4d6f36c61098de06-16-F', 'RR-02-a4e51784e182a017-16-F', 'RR-02-1ff6dce7b962afe3-16-F', 'RR-02-76e374dd28d81248-16-F', 'RR-02-3f7be56e3f9f187f-16-F', 'RR-02-ce0036d18be1cbb7-16-F', 'RR-02-d386eb98d569f14f-16-F', 'RR-02-19b55fc4074ea7a7-16-F', 'RR-02-993c9253f78c19bf-16-F', 'RR-02-4d30d6315e8e8d0e-16-F', 'RR-02-69bd577ca2f4d10e-16-F', 'RR-02-6be59bc48287ff0a-16-F', 'RR-02-a609d1ef085df966-16-F', 'RR-02-91911ecdb81cbdbc-16-F', 'RR-02-094f718f9fbf3087-16-F', 'RR-02-65c0adfb547ca69b-16-F', 'RR-02-8818dea65ff31cc5-16-F', 'RR-02-dc579a8e1f8ce4a0-16-F', 'RR-02-af8d6139571f9b75-16-F', 'RR-02-e810bb4aa2b8d414-16-F', 'RR-02-92ea06368d8a1330-16-F', 'RR-02-7a86b7668b9e0fc8-16-F', 'RR-02-ac6d35e02a432ef2-16-F', 'RR-02-3b31aefba79f28f5-16-F', 'RR-02-1e8e39bd86441957-16-F', 'RR-02-721183f2c545804f-16-F', 'RR-02-8b4324fcfd2b0746-16-F', 'RR-02-6fec6966f1f28405-16-F', 'RR-02-b618db3f535a34a5-16-F', 'RR-02-5bed8bcb898a8d24-16-F', 'RR-02-c89a1f53ed771996-16-F', 'RR-02-340f1860db055af1-16-F', 'RR-02-6f335b33a0a1be9d-16-F', 'RR-02-466420eb1f9fce92-16-F', 'RR-02-b052b51fefa9ce67-16-F', 'RR-02-8c1d80d26269af63-16-F', 'RR-02-6122c11077a9006c-16-F', 'RR-02-38bbc744ca2ac0b7-16-F', 'RR-02-f75e0e113579efa6-16-F', 'RR-02-47ecc93f52bc6e0b-16-F', 'RR-02-53f7dd11428601d8-16-F', 'RR-02-b87dd759d81bf89d-16-F', 'RR-02-3f6b8212b11c301c-16-F', 'RR-02-17f237b35b9c36c7-16-F', 'RR-02-ba05f56e7531b8ac-16-F', 'RR-02-c469571bf5dd3260-16-F', 'RR-02-757af66bae9301da-16-F', 'RR-02-0be9b0efc9e1b2d8-16-F', 'RR-02-3d3a456f054a2c75-16-F', 'RR-02-db0cdb3811bef3dd-16-F', 'RR-02-4567117770408ba1-16-F', 'RR-02-234e5a367f8d232a-16-F', 'RR-02-c3ec74eeb9cbf299-16-F', 'RR-02-80c1be875426e6d6-16-F', 'RR-02-fec67f6f4132cc0f-16-F', 'RR-02-14c9ef7ac8a2dc03-16-F', 'RR-02-2ab116f84d68e043-16-F', 'RR-02-469888aee72fbb8a-16-F', 'RR-02-c2c66cedf094533b-16-F', 'RR-02-cc9f5226009ba6e4-16-F', 'RR-02-289a74d34301b7af-16-F', 'RR-02-b62f37c6a404e83b-16-F', 'RR-02-e68f663510b3eadd-16-F', 'RR-02-bfbf4361736f996c-16-F', 'RR-02-d71acaf8e1fdd557-16-F', 'RR-02-1c5e2dd6014de992-16-F', 'RR-02-fe69aaf88d7c454a-16-F', 'RR-02-b8b5ef05dd6f8e12-16-F', 'RR-02-b4b4e7eff745e093-16-F', 'RR-02-7508e4dad8539dd3-16-F', 'RR-02-60c1aab00243be07-16-F', 'RR-02-052f180261baa23b-16-F', 'RR-02-2355de6b0a7f41ac-16-F', 'RR-02-73f4f0b40f22ed86-16-F', 'RR-02-42464e412b6a0366-16-F', 'RR-02-a27c101216726bd7-16-F', 'RR-02-4004b7ea811f0604-16-F', 'RR-02-2704be5eea79b7ae-16-F', 'RR-02-7e4de265f24b8a72-16-F', 'RR-02-d60e6efb94fb68b8-16-F', 'RR-02-40f28277bfc5ea8e-16-F', 'RR-02-4fe900904d0a9c85-16-F', 'RR-02-36bb1f194ab4866c-16-F', 'RR-02-cae9afc75fedcd22-16-F', 'RR-02-a15388df9112c5a6-16-F', 'RR-02-afeda5dd9e00d9da-16-F', 'RR-02-2fa9f69671f306d2-16-F', 'RR-02-5d459efbf36a1961-16-F', 'RR-02-0b06353b386ca535-16-F', 'RR-02-eddc4c82bc7f81b9-16-F', 'RR-02-d56f2788cec92217-16-F', 'RR-02-aa26c97ac100aa1a-16-F', 'RR-02-5549dac72492db6a-16-F', 'RR-02-03b5810baa310776-16-F', 'RR-02-3ba0719fa57afbea-16-F', 'RR-02-7c7ab7ac6b98a074-16-F', 'RR-02-b8b0271c31114118-16-F', 'RR-02-71834150ed737cbd-16-F', 'RR-02-6454a14b930bbb33-16-F', 'RR-02-b9ad3b92215e315a-16-F', 'RR-02-f2fb1c5efcb7140f-16-F', 'RR-02-bcd763d5f2e046be-16-F', 'RR-02-f59041f27577b37f-16-F', 'RR-02-860446453e33268d-16-F', 'RR-02-10c14f73d7c8d87e-16-F', 'RR-02-b76dd8e8d3e128df-16-F', 'RR-02-910859fe14f15abd-16-F', 'RR-02-b8ac6a83aeabdedd-16-F', 'RR-02-c8d7128110cc42d5-16-F', 'RR-02-ad1fcf135fe71ea3-16-F', 'RR-02-1c65c477fa187bde-16-F', 'RR-02-5329306871553552-16-F', 'RR-02-ab535f7af6f93fe1-16-F', 'RR-02-d3ff173f19d968c7-16-F', 'RR-02-e14ce618aee8279a-16-F', 'RR-02-085db4297707abff-16-F', 'RR-02-c5375175a751e58d-16-F', 'RR-02-031bd961221f7136-16-F', 'RR-02-40baa9f9769c49bd-16-F', 'RR-02-dd377c6c9f41d292-16-F'
]

ethylene_keys_archive = ['MNXR10316_MNXM977', 'MNXR14888_MNXM3270', 'MNXR16111_MNXM2550', 'MNXR16150_MNXM977', 'MNXR16310_MNXM2550', 'MNXR16310_MNXM3270', 'MNXR25094_MNXM2550', 'MNXR25094_MNXM977', 'MNXR77370_MNXM3270', 'MNXR8187_MNXM977', 'MNXR85539_MNXM977', 'MNXR87706_MNXM977', 'MNXR93892_MNXM977', 'MNXR9562_MNXM977']

caroten_keys = [
'RR-01-4ff35efaba69-12-F_12', 'RR-01-441413e93dd9-12-F_12', 'RR-01-e05fa6c58030-12-F_12',
'RR-01-e159cf6aa074-12-F_12', 'RR-01-1a456abbc0d8-12-F_12', 'RR-01-8910fd3fbc24-12-F_12',
'RR-01-1178f64b971e-12-F_12', 'RR-01-6d3500208261-12-F_12', 'RR-01-79a67c44bcbb-12-F_12',
'RR-01-44096cd01e46-12-F_12', 'RR-01-ab84cf3ac489-12-F_12', 'RR-01-693dcfc0cb83-12-F_12',
'RR-01-7c39c8a09b89-12-F_12', 'RR-01-1498c473c4cc-12-F_12', 'RR-01-6d2e028722f0-12-F_12',
'RR-01-35f0b9c12d05-12-F_12', 'RR-01-f7aa682a3ace-12-F_12', 'RR-01-91d8d1cbf1c1-12-F_12',
'RR-01-4daf335327a9-12-F_12', 'RR-01-ed22660e3711-12-F_12', 'RR-01-c9e94bc17bd5-12-F_12',
'RR-01-46b6c7e92d01-12-F_12', 'RR-01-9792e497c4e7-12-F_12', 'RR-01-a00d34ecee17-12-F_12',
'RR-01-38fbd0308afe-12-F_12', 'RR-01-e58379f37a1a-12-F_12', 'RR-01-20dadc0fd5e8-12-F_12',
'RR-01-5ef5450abd00-12-F_12', 'RR-01-20bda621c732-12-F_12', 'RR-01-c68691baf417-12-F_12',
'RR-01-0c7bceb70e55-12-F_12', 'RR-01-dd038e613058-12-F_12', 'RR-01-c687e462e254-12-F_12',
'RR-01-081e6a581b75-12-F_12', 'RR-01-0351542ed10f-12-F_12', 'RR-01-90bcfd187709-12-F_12',
'RR-01-2be959c35eb2-12-F_12', 'RR-01-11cb6d5fa275-12-F_12']

caroten_16_keys = [
    'RR-02-081e6a581b75d3b9-10-F', 'RR-02-2be959c35eb2c6ac-16-F', 'RR-02-11cb6d5fa275d4ce-10-F', 'RR-02-4ff35efaba6994fb-12-F', 'RR-02-441413e93dd96e37-12-F', 'RR-02-e05fa6c5803015b4-12-F', 'RR-02-e159cf6aa0744228-12-F', 'RR-02-1a456abbc0d8ce54-12-F', 'RR-02-8910fd3fbc24b470-12-F', 'RR-02-1178f64b971e36da-12-F', 'RR-02-6d35002082610a76-12-F', 'RR-02-79a67c44bcbbcf69-12-F', 'RR-02-44096cd01e46c9cf-12-F', 'RR-02-693dcfc0cb83c18e-12-F', 'RR-02-7c39c8a09b89f9f9-12-F', 'RR-02-1498c473c4cc62e9-12-F', 'RR-02-6d2e028722f03c3a-12-F', 'RR-02-35f0b9c12d0571b1-12-F', 'RR-02-f7aa682a3ace851c-12-F', 'RR-02-91d8d1cbf1c10bb6-12-F', 'RR-02-4daf335327a9d1b6-12-F', 'RR-02-ed22660e3711b360-12-F', 'RR-02-c9e94bc17bd5f6a3-12-F', 'RR-02-46b6c7e92d015b1b-12-F', 'RR-02-a00d34ecee17155c-12-F', 'RR-02-38fbd0308afee065-12-F', 'RR-02-e58379f37a1ae3ee-12-F', 'RR-02-5ef5450abd006d5a-12-F', 'RR-02-20bda621c732d5a9-12-F', 'RR-02-c68691baf417f77c-12-F', 'RR-02-0c7bceb70e55f042-12-F', 'RR-02-dd038e613058ed7e-12-F', 'RR-02-c687e462e2545683-12-F', 'RR-02-0351542ed10fdadf-12-F', 'RR-02-90bcfd187709eb1b-16-F'
]
caroten_keys_archive = ['MNXR15454_MNXM614', 'MNXR15602_MNXM557', 'MNXR16029_MNXM1048', 'MNXR17616_MNXM10398', 'MNXR17616_MNXM3549', 'MNXR26357_MNXM445', 'MNXR26357_MNXM557', 'MNXR26773_MNXM949', 'MNXR26773_MNXM557', 'MNXR26774_MNXM2099', 'MNXR26774_MNXM557', 'MNXR55435_MNXM139', 'MNXR60336_MNXM2000', 'MNXR70495_MNXM557', 'MNXR71200_MNXM1788', 'MNXR8603_MNXM1910', 'MNXR87812_MNXM557']


def get_rules_and_score(full_rules_forward_H = None,
                        full_rules_retro_H = None,
                        full_rules_forward_no_H = None,
                        full_rules_retro_no_H = None,
                        add_Hs = True,
                        retro = True,
                        diameters = [16],
                        small = False,
                        c_name = None,
                        filtering_EC = None):
    """
    full_rules: fondation data
        - H or no H
        - forward or reverse
        Those rules will be calculated in calculate_rules_similarity for those 4 conditions.
        Expert users can choose to calculate their own rules under the required format
        and specify only those they intedn to use.
    diameter: a list of diameters to use (it's a rule attribute)
    small: current optional argument for running fast tests.
    """
    if add_Hs:
        if retro:
            unfiltered_rules = full_rules_retro_H
            scoring = "BiologicalFullScoringRetroH"
        else:
            unfiltered_rules = full_rules_forward_H
            scoring = "BiologicalFullScoringFwdH"
    else:
        if retro:
            unfiltered_rules = full_rules_retro_no_H
            scoring = "BiologicalFullScoringRetroNoH"
        else:
            unfiltered_rules = full_rules_forward_no_H
            scoring = "BiologicalFullScoringFwdNoH"

    if small:
        merged_rules = full_rules_retro_H
        merged_rules.update(full_rules_forward_H)
        if c_name == "TPA":
            rules = dict((rid, merged_rules[rid]) for rid in TPA_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "deoxiviolacein" or c_name == "violacein":
            rules = dict((rid, merged_rules[rid]) for rid in violacein_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "pinocembrin":
            rules = dict((rid, merged_rules[rid]) for rid in pinocembrin_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "naringenin":
            rules = dict((rid, merged_rules[rid]) for rid in naringenin_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "ethylene":
            rules = dict((rid, merged_rules[rid]) for rid in ethylene_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "caroten":
            rules = dict((rid, merged_rules[rid]) for rid in caroten_16_keys)
            scoring = "BiologicalFullScoringH"
        elif c_name == "styrene":
            rules = dict((rid, merged_rules[rid]) for rid in styrene_16_keys)
            scoring = "BiologicalFullScoringH"
        else:
            if biosensor:
                rules = unfiltered_rules
            else:
                raise NotImplementedError
        return(rules, scoring)
    else:
        # Filtering by diameters
        if diameters is None:
            filtered_rules = unfiltered_rules
        else:
            filtered_rules = {}
            for rid, rule_characteritics in unfiltered_rules.items():
                if int(rule_characteritics["Diameter"]) in diameters:
                    filtered_rules[rid] = rule_characteritics
        # Filtering by EC:
        if filtering_EC is None:
            filtered_rules_EC = filtered_rules
        else:
            filtered_rules_EC = {}
            for rid, rule_characteritics in filtered_rules.items():
                for ec_reac in rule_characteritics["EC_number"]:
                    for ec_query in filtering_EC:
                        ec_query_split = ec_query.split('.')
                        ec_correct_depth = get_EC_from_depth(ec_reac, len(ec_query_split))
                        if ec_correct_depth == ec_query:
                            filtered_rules_EC[rid] = rule_characteritics

        return(filtered_rules_EC, scoring)

def get_EC_from_depth(EC, EC_depth = 4):
    """
    Obtains an EC number at the requested depth from a full EC number.
    Returns None is the requested depth is too long (ie: 4.1.-- would return None at depth 3 or 4)
    """
    EC_split = EC.split('.')
    if len(EC_split) < EC_depth:
        return(None)
    else:
        EC_cut = ".".join([EC_split[i] for i in range(EC_depth)])
        return(EC_cut)
