"""
Defines organisms as chemical_compounds_state objects.
Unpickled after calculation when setting up RP3.
"""

# General utilities
import logging
import pickle
import os
import csv
import sys

from config import *

# RP3 specific objects
from compound import Compound, unpickle
from chemical_compounds_state import ChemicalCompoundState
from rdkit.Chem import AllChem
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError
from utilities.reactor.cli import worker_match, worker_fire, kill, RuleConversionError


class NotReady(Exception):
    """Raised when organisms or rules have not been caculated in advance"""

    def __init__(self, msg = "Not Ready. Need to run set-up scripts"):
        self._msg = msg

    def __str__(self):
        return self._msg

def import_organism_from_csv(csv_file, add_Hs = True):
    with open(csv_file) as csv_handle:
        dict_reader = csv.DictReader(csv_handle, delimiter = ',')
        compound_list = []
        for row in dict_reader:
            name = row["name"]
            inchi = row["inchi"]
            if inchi is None or inchi == "None" or inchi =="":
                pass
            else:
                try:
                    if name.startswith("InChI"):
                        compound = Compound(InChI = inchi, heavy_standardisation = True, force_add_H = add_Hs)
                    else:
                        compound = Compound(InChI = inchi, name = name, heavy_standardisation = True, force_add_H = add_Hs)
                    if not compound.in_list(compound_list, main_layer = False):
                        compound_list.append(compound)
                except ChemConversionError as e:
                    logging.error("For compound {} with inchi {}: error ChemConversionError".format(name, inchi))
    organism = ChemicalCompoundState(compound_list, main_layer = False)
    # organism.set_main_layer(True)
    return(organism)


organisms_data_path = "{}/organisms".format(data_path)
if not os.path.exists(organisms_data_path):
    os.mkdir(organisms_data_path)

if not os.path.exists(organisms_data_path + '/state_iML1515_chassis_H.pkl'):
    logging.error("Please run calculate_organisms script")
    raise NotReady


Test_organism_H = unpickle(file_name = "{}".format('Test_organism_H'), type = 'state', folder_address = organisms_data_path)
ecoli_chassis_H = unpickle(file_name = "{}".format('iML1515_chassis_H'), type = 'state', folder_address = organisms_data_path)
detectable_cmpds_H = unpickle(file_name = "{}".format('detectable_cmpds_H'), type = 'state', folder_address = organisms_data_path)
core_ecoli_H = unpickle(file_name = "{}".format('core_ecoli_H'), type = 'state', folder_address = organisms_data_path)
bsubtilis_H = unpickle(file_name = "{}".format('bsubtilis_H'), type = 'state', folder_address = organisms_data_path)
iJO1366_chassis_H = unpickle(file_name = "{}".format('iJO1366_chassis_H'), type = 'state', folder_address = organisms_data_path)


Test_organism_noH = unpickle(file_name = "{}".format('Test_organism_noH'), type = 'state', folder_address = organisms_data_path)
ecoli_chassis_noH = unpickle(file_name = "{}".format('iML1515_chassis_noH'), type = 'state', folder_address = organisms_data_path)
detectable_cmpds_noH = unpickle(file_name = "{}".format('detectable_cmpds_noH'), type = 'state', folder_address = organisms_data_path)
core_ecoli_noH = unpickle(file_name = "{}".format('core_ecoli_noH'), type = 'state', folder_address = organisms_data_path)
bsubtilis_noH = unpickle(file_name = "{}".format('bsubtilis_noH'), type = 'state', folder_address = organisms_data_path)
iJO1366_chassis_noH = unpickle(file_name = "{}".format('iJO1366_chassis_noH'), type = 'state', folder_address = organisms_data_path)
