"""
This module loads calculates organisms
- standardises compounds within the organism
- saves them as pickles that can be laoded by RP3
"""

# General utilities
import logging
import pickle
import os
import csv
import sys
import argparse

from config import *

# RP3 specific objects
from compound import Compound, unpickle
from chemical_compounds_state import ChemicalCompoundState
from rdkit.Chem import AllChem
from utilities.reactor.Utils import standardize_chemical, standardize_results, handle_results, ChemConversionError
from utilities.reactor.cli import worker_match, worker_fire, kill, RuleConversionError

def __run__():
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
        return(organism)

    # Calculate with H
    directory_for_saving = "{}/{}".format(os.getcwd(), os.path.dirname(__file__))
    compound_1 = Compound("[H+]", name = "1", heavy_standardisation = True, force_add_H = True)
    compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", force_add_H = True, name = '6', heavy_standardisation = True)
    compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]", name = '3459', heavy_standardisation = True, force_add_H = True)
    Test_organism = ChemicalCompoundState(state_name = "Test", compound_list = [compound_1, compound_6, compound_3459])
    detectable_cmpds = import_organism_from_csv('{}/data/sinks/detectable_metabolites_uncommented.csv'.format(directory_for_saving), add_Hs = True)
    iML1515_chassis = import_organism_from_csv("{}/data/sinks/ecoli_iML1515_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = True)
    core_ecoli = import_organism_from_csv("{}/data/sinks/ecoli_core_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = True)
    iJO1366_chassis = import_organism_from_csv("{}/data/sinks/ecoli_iJO1366_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = True)
    bsubtilis = import_organism_from_csv("{}/data/sinks/bsubtilis_iYO844_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = True)
    Test_organism.save(file_name = 'Test_organism_H', folder_address = organisms_data_path)
    detectable_cmpds.save(file_name = 'detectable_cmpds_H', folder_address = organisms_data_path)
    iML1515_chassis.save(file_name = 'iML1515_chassis_H', folder_address = organisms_data_path)
    core_ecoli.save(file_name = "{}".format('core_ecoli_H'), folder_address = organisms_data_path)
    iJO1366_chassis.save(file_name = "{}".format('iJO1366_chassis_H'), folder_address = organisms_data_path)
    bsubtilis.save(file_name = "{}".format('bsubtilis_H'), folder_address = organisms_data_path)
    # Calculate without H
    compound_1 = Compound("[H+]", name = "1", heavy_standardisation = True, force_add_H = False)
    compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", force_add_H = False, name = '6', heavy_standardisation = True)
    compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]", name = '3459', heavy_standardisation = True, force_add_H = False)
    Test_organism = ChemicalCompoundState(state_name = "Test", compound_list = [compound_1, compound_6, compound_3459])
    detectable_cmpds = import_organism_from_csv('{}/data/sinks/detectable_metabolites_uncommented.csv'.format(directory_for_saving), add_Hs = False)
    iML1515_chassis = import_organism_from_csv("{}/data/sinks/ecoli_iML1515_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = False)
    core_ecoli = import_organism_from_csv("{}/data/sinks/ecoli_core_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = False)
    iJO1366_chassis = import_organism_from_csv("{}/data/sinks/ecoli_iJO1366_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = False)
    bsubtilis = import_organism_from_csv("{}/data/sinks/bsubtilis_iYO844_sink_reduced_rp_ready.csv".format(directory_for_saving), add_Hs = False)
    Test_organism.save(file_name = 'Test_organism_noH', folder_address = organisms_data_path)
    detectable_cmpds.save(file_name = 'detectable_cmpds_noH', folder_address = organisms_data_path)
    iML1515_chassis.save(file_name = 'iML1515_chassis_noH', folder_address = organisms_data_path)
    core_ecoli.save(file_name = "{}".format('core_ecoli_noH'), folder_address = organisms_data_path)
    iJO1366_chassis.save(file_name = "{}".format('iJO1366_chassis_noH'), folder_address = organisms_data_path)
    bsubtilis.save(file_name = "{}".format('bsubtilis_noH'), folder_address = organisms_data_path)


if __name__ == "__main__":
    d = "Formatting organisms in a RP3 compatible format"
    parser = argparse.ArgumentParser(description=d)
    parser.add_argument("--terminal", help="Default logger is logs_organisms_set_up, switch to terminal if specified",
                        action='store_true', default=False)
    args = parser.parse_args()
    organisms_data_path = "{}/organisms".format(data_path)
    if not os.path.exists(organisms_data_path):
        os.mkdir(organisms_data_path)
    if args.terminal is True:
        logging.basicConfig(
                stream = sys.stderr,
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
    else:
        logging.basicConfig(
                stream = open("{}/{}.log".format(organisms_data_path, "logs_organisms_set_up"), "w"),
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
        print("By default, logs are saved in {}/logs_organisms_set_up.log. Please use --terminal to redirect to sys.stderr".format(organisms_data_path))
    __run__()
