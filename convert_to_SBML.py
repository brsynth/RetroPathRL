"""
Converts pathways under json format to SBML format
"""

# General utilities
import sys
import logging
import csv
import copy
import json
import pickle
import libsbml
from hashlib import md5
import os
import argparse

# RP3 specific objects
from compound import Compound
from move import Move

def _nameToSbmlId(name):
    IdStream = []
    count = 0
    end = len(name)
    if '0' <= name[count] and name[count] <= '9':
        IdStream.append('_')
    for count in range(0, end):
        if (('0' <= name[count] and name[count] <= '9') or
                ('a' <= name[count] and name[count] <= 'z') or
                ('A' <= name[count] and name[count] <= 'Z')):
            IdStream.append(name[count])
        else:
            IdStream.append('_')
    Id = ''.join(IdStream)
    if Id[len(Id) - 1] != '_':
        return Id
    return Id[:-1]

def add_specy(sbml_model,
             chemId = 'Id_cmpound',
             smiles = "smilescomppoun",
             inchi = "inchicompounds",
             inchiKey = "inchiKeycomppoun",
             name = "compounds_name",
             in_sink = False):

    spe = sbml_model.createSpecies()
    spe.setCompartment("cytoplasm")
    spe.setHasOnlySubstanceUnits(False)
    spe.setBoundaryCondition(False)
    spe.setConstant(False)
    spe.setInitialConcentration(1.0)
    clean_id = str(chemId)+'__64__'+str("cytoplasm")
    clean_id = clean_id.replace('-', '_')  # No - in name
    metaid = _nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())
    spe.setMetaId(metaid)
    spe.setName(name)
    if in_sink:
        annotation = '''<annotation>
        <rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'''
        annotation += '''
        <rdf:RP3 rdf:about="#'''+str(metaid or '')+'''">
        <RP3:RP3 xmlns:RP3="https://github.com/brsynth">
        <RP3:smiles>'''+str(smiles or '')+'''</RP3:smiles>
        <RP3:inchi>'''+str(inchi or '')+'''</RP3:inchi>
        <RP3:inchikey>'''+str(inchiKey or '')+'''</RP3:inchikey>
        <RP3:in_sink>'''+ str(True)+'''</RP3:in_sink>
        </RP3:RP3>
        </rdf:RP3>'''
        annotation += '''
        </rdf:RDF>
        </annotation>'''
    else:
        annotation = '''<annotation>
        <rdf:RDF
        xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'''
        annotation += '''
        <rdf:RP3 rdf:about="#'''+str(metaid or '')+'''">
        <RP3:RP3 xmlns:RP3="https://github.com/brsynth">
        <RP3:smiles>'''+str(smiles or '')+'''</RP3:smiles>
        <RP3:inchi>'''+str(inchi or '')+'''</RP3:inchi>
        <RP3:inchikey>'''+str(inchiKey or '')+'''</RP3:inchikey>
        </RP3:RP3>
        </rdf:RP3>'''
        annotation += '''
        </rdf:RDF>
        </annotation>'''
    spe.setAnnotation(annotation)
    return(sbml_model)

def add_reaction(sbml_model,
                 reacId = 'Id_reac',
                 ec = "Test_ec",
                 rule_id = "rule_id",
                 biological_score = "biological_score",
                 chemical_score = "chemical_score",
                 reactant_stoechio = {},
                 product = "product_name",
                 reaction_smiles = "reaction_smiles",
                 diameter = "diameter"):
    reac = sbml_model.createReaction()

    reac_fbc = reac.getPlugin('fbc')
    reac_fbc.setUpperFluxBound('B_999999')
    reac_fbc.setLowerFluxBound('B_0')
    #reactions
    reac.setId(reacId)
    reac.setSBOTerm(185)
    reac.setReversible(True)
    reac.setFast(False)
    metaid = _nameToSbmlId(md5(str(reacId).encode('utf-8')).hexdigest())
    reac.setMetaId(metaid)
    #reactants_dict
    for reactant in reactant_stoechio.keys():
        chemId = reactant
        spe = reac.createReactant()
        clean_id = str(chemId)+'__64__'+str("cytoplasm")
        clean_id = clean_id.replace('-', '_')  # No - in name
        spe.setSpecies(clean_id)
        spe.setConstant(True)
        try:
            stoechio = reactant_stoechio[reactant]
        except KeyError:
            stoechio = 1
        spe.setStoichiometry(stoechio)
    #products_dict
    if not product is None:
        pro = reac.createProduct()
        clean_id = str(product)+'__64__'+str("cytoplasm")
        clean_id = clean_id.replace('-', '_')  # No - in name
        pro.setSpecies(clean_id)
        pro.setConstant(True)
        pro.setStoichiometry(1)
    #annotation
    annotation = '''<annotation>
    <rdf:RDF
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#">'''

    annotation += '''
    <rdf:RP3 rdf:about="#'''+str(metaid or '')+'''">
      <RP3:RP3 xmlns:RP3="https://github.com/brsynth">
        <RP3:smiles>'''+str(reaction_smiles or '')+'''</RP3:smiles>
        <RP3:rule_id>'''+str(rule_id or '')+'''</RP3:rule_id>
        <RP3:EC>'''+str(ec)+'''</RP3:EC>
        <RP3:biological_score value="'''+str(biological_score or '')+'''" />
        <RP3:chemical_score value="'''+str(chemical_score or '')+'''" />
        <RP3:diameter value="'''+str(diameter or '')+'''" />
      </RP3:RP3>
    </rdf:RP3>
    </rdf:RDF>
    </annotation>'''
    reac.setAnnotation(annotation)
    return(sbml_model)


def convert_json_to_SBML(json_file, modelID = "test", folder_to_save = 'temp'):
    # Set up the empty model
    smbl_namespace = libsbml.SBMLNamespaces(3,1)
    smbl_namespace.addPkgNamespace('fbc',2)
    smbl_namespace.addPkgNamespace('groups',2)
    document = libsbml.SBMLDocument(smbl_namespace)
    sbml_model = document.createModel()
    sbml_model.getPlugin('fbc')
    sbml_model.getPlugin('groups')
    sbml_model.setId(modelID)
    sbml_model.setName(modelID)
    sbml_model.setTimeUnits('second')
    sbml_model.setExtentUnits('mole')
    sbml_model.setSubstanceUnits('mole')
    # Could implement units, currently removed from the model
    # Should have it in a seperate function
    compartment = sbml_model.createCompartment()
    compartment.setId("cytoplasm")
    target_node = None
    for node in json_file["elements"]["nodes"]:
        if node["data"]["type"] == "compound":
            sbml_model = add_specy(sbml_model,
                      chemId = node["data"]["id"],
                      smiles = node["data"]["SMILES"],
                      inchi = node["data"]["InChI"],
                      inchiKey = node["data"]["id"],
                      name = ",".join(node["data"]["Names"]),
                      in_sink = node["data"]["inSink"] == 1)
        if node["data"]["isSource"] == 1:
            logging.info("Target node is {}".format(node["data"]["id"]))
            target_node = node
    for element in sbml_model.getListOfSpecies():
        logging.debug(element)
    for node in json_file["elements"]["nodes"]:
        if node["data"]["type"] == "reaction":
            try:
                reactant_stoechio = node["data"]["Stoechiometry"]
            except KeyError:
                reactant_stoechio = {}
            sbml_model = add_reaction(sbml_model,
                                     reacId = node["data"]["id"],
                                     ec = ','.join(node["data"]["EC number"]),
                                     rule_id = ','.join(node["data"]["Rule ID"]),
                                     biological_score = node["data"]["Score"],
                                     chemical_score = node["data"]["ChemicalScore"],
                                     reactant_stoechio = reactant_stoechio,
                                     product = node["data"]["id"].split("-RR")[0],
                                     reaction_smiles = node["data"]["Reaction SMILES"],
                                     diameter = node["data"]["Diameter"])
    sbml_model = add_reaction(sbml_model,
                             reacId = "production",
                             ec = 'NA',
                             rule_id = 'NA',
                             biological_score = 'NA',
                             chemical_score = 'NA',
                             reactant_stoechio = {target_node["data"]["id"]: 1},
                             product = None,
                             reaction_smiles = 'NA',
                             diameter = 'NA')

    document.setModel(sbml_model)
    libsbml.writeSBMLToFile(document,'{}/{}.xml'.format(folder_to_save, modelID))
    pass


def __cli():
    def define_folder_to_save(folder):
        if folder is None:
            folder_to_save = os.path.join('debugging_results', args.c_name)
        else:
            folder_to_save = folder
        if not os.path.exists(folder_to_save):
            os.makedirs(folder_to_save, exist_ok=True)
        return folder_to_save
    d = "Command line interface to convert json files to SBML files"
    parser = argparse.ArgumentParser(description=d)
    # Logs and saving information
    """Command line interface to convert json files to SBML files"""
    parser.add_argument("--verbose", help="Default logger is INFO, switch to DEBUG is specified",
                        dest='verbose', action='store_true', default=False)
    parser.add_argument("--log_file", help="Default logger is stderr, switch to log_file if specified",
                        default=None)
    parser.add_argument("--folder_to_save",
                        help="Folder to store results. Default: temp",
                        default="temp")
    parser.add_argument("--json_convert",
                        help="File to convert",
                        default="deoxi_07_no_H/deoxiviolacein_iteration_85.json")
    parser.add_argument("--file_name", help = 'File name if name changes.', default = None)
    args = parser.parse_args()
    # Setting up the logs
    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    if args.log_file is None:
        logging.basicConfig(stream=sys.stderr,
                            level=logging_level,
                            datefmt='%d/%m/%Y %H:%M:%S',
                            format='%(asctime)s -- %(levelname)s -- %(message)s')
    else:
        if not "log" in args.log_file:
            log_file = "log_" + args.log_file
        else:
            log_file = args.log_file
        log_writer = open("{}/{}".format(folder_to_save, log_file), "w")
        logging.basicConfig(stream=log_writer,
                            level=logging_level,
                            datefmt='%d/%m/%Y %H:%M:%S',
                            format='%(asctime)s -- %(levelname)s -- %(message)s')

    folder_to_save = define_folder_to_save(args.folder_to_save)
    # Choosing file
    if args.file_name is None:
        model_ID = args.json_convert.split("/")[-1].split(".json")[0]
    else:
        model_ID = args.file_name
    pathway_to_test = json.load(open(args.json_convert, "r"))
    convert_json_to_SBML(pathway_to_test, model_ID, folder_to_save = folder_to_save)

if __name__ == "__main__":
    __cli()
