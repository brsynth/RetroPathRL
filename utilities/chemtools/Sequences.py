#!/usr/bin/env python
"""Sequences of filters to be used for standardization."""


from utilities.chemtools.Filters import Filters
from rdkit.Chem import Cleanup, SanitizeMol, SanitizeFlags
from rdkit.Chem.AllChem import AssignStereochemistry


def sequence_rr_legacy(mol):
    """Sequence of filters applied for the first version of RetroRules
    """
    F = Filters()
    Cleanup(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
    mol = F.remove_isotope(mol)
    mol = F.neutralise_charge(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    mol = F.keep_biggest(mol)
    mol = F.add_hydrogen(mol, addCoords=True)
    mol = F.kekulize(mol)
    return mol


def sequence_tunable(
        mol,
        OP_REMOVE_ISOTOPE=True, OP_NEUTRALISE_CHARGE=True,
        OP_REMOVE_STEREO=False, OP_COMMUTE_INCHI=False,
        OP_KEEP_BIGGEST=True, OP_ADD_HYDROGEN=True,
        OP_KEKULIZE=True, OP_NEUTRALISE_CHARGE_LATE=True
    ):
    """Tunable sequence of filters for standardization.
    
    Operations will made in the following order:
     1 RDKit Cleanup      -- always
     2 RDKIT SanitizeMol  -- always
     3 Remove isotope     -- optional (default: True)
     4 Neutralise charges -- optional (default: True)
     5 RDKit SanitizeMol  -- if 4 or 5
     6 Remove stereo      -- optional (default: False)
     7 Commute Inchi      -- if 6 or optional (default: False)
     8 Keep biggest       -- optional (default: True)
     9 RDKit SanitizeMol  -- if any (6, 7, 8)
    10 Add hydrogens      -- optional (default: True)
    11 Kekulize           -- optional (default: True)
    """
    F = Filters()
    # Always perform the basics..
    Cleanup(mol)
    SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
    # 
    if OP_REMOVE_ISOTOPE:
        mol = F.remove_isotope(mol)
    if OP_NEUTRALISE_CHARGE:
        mol = F.neutralise_charge(mol)
    if any([OP_REMOVE_ISOTOPE, OP_REMOVE_ISOTOPE]):
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    # 
    if OP_REMOVE_STEREO:
        mol = F.remove_stereo(mol)
        OP_COMMUTE_INCHI = True
    if OP_COMMUTE_INCHI:
        mol = F.commute_inchi(mol)
    if OP_KEEP_BIGGEST:
        mol = F.keep_biggest(mol)
    if any([OP_REMOVE_STEREO, OP_COMMUTE_INCHI, OP_KEEP_BIGGEST]):
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    #
    if OP_NEUTRALISE_CHARGE_LATE:
        mol = F.neutralise_charge(mol)
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
    #
    if OP_ADD_HYDROGEN:
        mol = F.add_hydrogen(mol, addCoords=True)
    if OP_KEKULIZE:
        mol = F.kekulize(mol)
    #
    return mol
