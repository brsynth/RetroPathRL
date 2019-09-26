#!/usr/bin/env python
"""
Starting a new toolbox to handle chemical compounds
"""

from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs


def convert_depiction(idepic, itype='smiles', otype={'inchikey'}):
    """Convert chemical depiction to others type of depictions
    
    :param  idepic: string depiction to be converted, str
    :param   itype: type of depiction provided as input, str
    :param   otype: types of depiction to be generated, {"", "", ..}
    :return odepic: generated depictions, {"otype1": "odepic1", ..}
    
    Usage example:
    - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    """
    # Import (if needed)
    if itype == 'smiles':
        rdmol = MolFromSmiles(idepic, sanitize=True)
    elif itype == 'inchi':
        rdmol = MolFromInchi(idepic, sanitize=True)
    else:
        raise NotImplementedError('"{}" is not a valid input type'.format(itype))
    if rdmol is None:  # Check imprt
        raise Exception('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
    
    # Export
    odepic = dict()
    for item in otype:
        if item == 'smiles':
            odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
        elif item == 'inchi':
            odepic[item] = MolToInchi(rdmol)
        elif item == 'inchikey':
            odepic[item] = MolToInchiKey(rdmol)
        else:
            raise NotImplementedError('"{}" is not a valid output type'.format(otype))

    return odepic
