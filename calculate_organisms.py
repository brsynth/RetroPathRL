"""
This module loads calculates organisms
- standardises compounds within the organism
- saves them as pickles that can be laoded by RP3
"""

# General utilities
import logging
import os
import csv
import sys
import argparse

from config import DATA_PATH

# RP3 specific objects
from compound import Compound
from chemical_compounds_state import ChemicalCompoundState
from utilities.reactor.Utils import ChemConversionError


def __run__():
    def import_organism_from_csv(csv_file, add_Hs = True):
        with open(csv_file) as csv_handle:
            dict_reader = csv.DictReader(csv_handle, delimiter = ',')
            compound_list = []
            for row in dict_reader:
                name = row["name"]
                inchi = row["inchi"]
                if inchi is None or inchi == "None" or inchi == "":
                    pass
                else:
                    try:
                        if name.startswith("InChI"):
                            compound = Compound(InChI = inchi, heavy_standardisation = True, force_add_H = add_Hs)
                        else:
                            compound = Compound(InChI = inchi, name = name, heavy_standardisation = True, force_add_H = add_Hs)
                        if not compound.in_list(compound_list, main_layer = False):
                            compound_list.append(compound)
                    except ChemConversionError:
                        logging.error("For compound {} with inchi {}: error ChemConversionError".format(name, inchi))
        organism = ChemicalCompoundState(compound_list, main_layer = False)
        return organism

    # Calculate with H ========================================================
    logging.info("Calculating organisms with H...")

    # Test organism
    compound_1 = Compound("[H+]", name="1", heavy_standardisation=True, force_add_H=True)
    compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", force_add_H=True, name='6', heavy_standardisation=True)
    compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]", name='3459', heavy_standardisation=True, force_add_H=True)
    test_organism = ChemicalCompoundState(state_name="Test", compound_list=[compound_1, compound_6, compound_3459])

    # Load real organisms
    detectable_cmpds = import_organism_from_csv(f"{SINK_DATA_PATH}/detectable_metabolites_uncommented.csv", add_Hs=True)
    iML1515_chassis = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_iML1515_sink_reduced_rp_ready.csv", add_Hs=True)
    core_ecoli = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_core_sink_reduced_rp_ready.csv", add_Hs=True)
    iJO1366_chassis = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_iJO1366_sink_reduced_rp_ready.csv", add_Hs=True)
    bsubtilis = import_organism_from_csv(f"{SINK_DATA_PATH}/bsubtilis_iYO844_sink_reduced_rp_ready.csv", add_Hs=True)

    # Save organisms
    test_organism.save(file_name='Test_organism_H', folder_address=ORGANISMS_DATA_PATH)
    detectable_cmpds.save(file_name='detectable_cmpds_H', folder_address=ORGANISMS_DATA_PATH)
    iML1515_chassis.save(file_name='iML1515_chassis_H', folder_address=ORGANISMS_DATA_PATH)
    core_ecoli.save(file_name='core_ecoli_H', folder_address=ORGANISMS_DATA_PATH)
    iJO1366_chassis.save(file_name='iJO1366_chassis_H', folder_address=ORGANISMS_DATA_PATH)
    bsubtilis.save(file_name='bsubtilis_H', folder_address=ORGANISMS_DATA_PATH)

    # Calculate without H =====================================================
    logging.info("Calculating organisms without H...")

    # Test organism
    compound_1 = Compound("[H+]", name="1", heavy_standardisation=True, force_add_H=False)
    compound_6 = Compound("[H][N]=[C]([O][H])[C]1=[C]([H])[N]([C]2([H])[O][C]([H])([C]([H])([H])[O][P](=[O])([O][H])[O][P](=[O])([O][H])[O][C]([H])([H])[C]3([H])[O][C]([H])([n]4[c]([H])[n][c]5[c]([N]([H])[H])[n][c]([H])[n][c]54)[C]([H])([O][P](=[O])([O][H])[O][H])[C]3([H])[O][H])[C]([H])([O][H])[C]2([H])[O][H])[C]([H])=[C]([H])[C]1([H])[H]", force_add_H=False, name='6', heavy_standardisation=True)
    compound_3459 = Compound("[H][O][C](=[O])[C](=[O])[C]([H])([H])[C]([H])([O][H])[C]([H])([O][H])[C]([H])([H])[H]", name='3459', heavy_standardisation=True, force_add_H=False)
    test_organism = ChemicalCompoundState(state_name="Test", compound_list=[compound_1, compound_6, compound_3459])

    # Load real organisms
    detectable_cmpds = import_organism_from_csv(f"{SINK_DATA_PATH}/detectable_metabolites_uncommented.csv", add_Hs=True)
    iML1515_chassis = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_iML1515_sink_reduced_rp_ready.csv", add_Hs=False)
    core_ecoli = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_core_sink_reduced_rp_ready.csv", add_Hs=False)
    iJO1366_chassis = import_organism_from_csv(f"{SINK_DATA_PATH}/ecoli_iJO1366_sink_reduced_rp_ready.csv", add_Hs=False)
    bsubtilis = import_organism_from_csv(f"{SINK_DATA_PATH}/bsubtilis_iYO844_sink_reduced_rp_ready.csv", add_Hs=False)

    # Save organisms
    test_organism.save(file_name='Test_organism_noH', folder_address=ORGANISMS_DATA_PATH)
    detectable_cmpds.save(file_name='detectable_cmpds_noH', folder_address=ORGANISMS_DATA_PATH)
    iML1515_chassis.save(file_name='iML1515_chassis_noH', folder_address=ORGANISMS_DATA_PATH)
    core_ecoli.save(file_name='core_ecoli_noH', folder_address=ORGANISMS_DATA_PATH)
    iJO1366_chassis.save(file_name='iJO1366_chassis_noH', folder_address=ORGANISMS_DATA_PATH)
    bsubtilis.save(file_name='bsubtilis_noH', folder_address=ORGANISMS_DATA_PATH)

    return 0


if __name__ == "__main__":
    d = "Formatting organisms in a RP3 compatible format"
    parser = argparse.ArgumentParser(description=d)
    parser.add_argument("--terminal", help="Default logger is logs_organisms_set_up, switch to terminal if specified",
                        action='store_true', default=False)
    args = parser.parse_args()

    # Sink data path
    global SINK_DATA_PATH
    SINK_DATA_PATH = f"{DATA_PATH}/sinks"
    assert os.path.exists(SINK_DATA_PATH), f"Sink data path {SINK_DATA_PATH} does not exist"

    # Organisms data path
    global ORGANISMS_DATA_PATH
    ORGANISMS_DATA_PATH = f"{DATA_PATH}/organisms"
    if not os.path.exists(ORGANISMS_DATA_PATH):
        os.mkdir(ORGANISMS_DATA_PATH)

    if args.terminal is True:
        logging.basicConfig(
                stream=sys.stderr,
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
    else:
        logging.basicConfig(
                stream=open("{}/{}.log".format(ORGANISMS_DATA_PATH, "logs_organisms_set_up"), "w"),
                level=logging.INFO,
                datefmt='%d/%m/%Y %H:%M:%S',
                format='%(asctime)s -- %(levelname)s -- %(message)s'
                )
        print(f"By default, logs are saved in {ORGANISMS_DATA_PATH}/logs_organisms_set_up.log. Please use --terminal to redirect to sys.stderr")
    __run__()
