#!/usr/bin/env python
"""
Standardize chemicals

This is basically a rework of the standardizer.py written by Baudoin Delepine
at INRA.

@author: Baudoin Delepine, 2016-2017
@author: Thomas Duigou, 2018-2019
"""

from utilities.chemtools import Sequences
from utilities.chemtools.Filters import Filters
from rdkit.Chem import SanitizeMol, SanitizeFlags
from rdkit.Chem.AllChem import AssignStereochemistry

class Standardizer(object):
    """Handle standardization of compound(s) through user-defined "filters".
    """

    def __call__(self, mol):
            """Calling the Standardizer class like a function is the same
            as calling its "compute" method.
            
            Form:
                https://github.com/mcs07/MolVS/blob/master/molvs/standardize.py
            """
            return self.compute(mol)

    def __init__(self, sequence_fun=None, params=None):
        """Set up parameters for the standardization
        
        :param rdmol: an RDKit Mol object
        """
        # Function to be used for standardizing compounds
        # Add you own function as method class
        if sequence_fun is None:
            self.sequence_fun = self.sequence_minimal
        elif callable(sequence_fun):  # Guess: fun_filters is the function itself
            self.sequence_fun = sequence_fun
        elif type(sequence_fun) == str:
            self.sequence_fun = getattr(Sequences, sequence_fun)  # Guess: sequence_fun is the name of the function
        # Arguments to be passed to any custom standardization function
        self._params = params if params else None

    def sequence_minimal(self, mol):
        """Minimal standardization."""
        SanitizeMol(mol, sanitizeOps=SanitizeFlags.SANITIZE_ALL, catchErrors=False)
        AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)  # Fix bug TD201904.01
        return mol
        
    def compute(self, mol):
        """Do the job."""
        if self._params is None:
            return self.sequence_fun(mol)
        else:
            return self.sequence_fun(mol, **self._params)
