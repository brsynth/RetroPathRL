#!/usr/bin/env python
"""
Set of filters to be used for chemical standardisation

@author: Baudoin Delepine, 2016-2017
@author: Thomas Duigou, 2018-2019
"""

from copy import deepcopy
from rdkit.Chem import AddHs, GetMolFrags, Kekulize, MolToInchi, MolFromInchi, MolFromSmarts, MolFromSmiles, RemoveStereochemistry, MolToSmiles, RemoveHs
from rdkit.Chem.AllChem import Compute2DCoords, ReplaceSubstructs
from rdkit.Chem.Descriptors import MolWt


class Filters(object):
    """Set of filters to be used for chemical standardization.
    """
    
    @classmethod
    def _copy_properties(cls, mol_from, mol_to):
        """Copy properties from a RDKit compound to another one.
        
        :param  mol_from: RDKit Mol source object
        :param  mol_to:   RDKit Mol target object

        Warning: aside from chemical's name, all private properties are lost.
        """
        # NB: name is stored in its default location which is "_Name" and
        # is a private propertie.
        property_list = mol_from.GetPropNames(includePrivate=False)
        if mol_from.HasProp('_Name'):  # TD: If _Name is set always save name
            property_list.append("_Name")
        for property_name in property_list:
            mol_to.SetProp(property_name, mol_from.GetProp(property_name))
    
    @classmethod
    def keep_biggest(cls, mol_in):
        """Strip small fragments from compound.

        Returns a new compound where only the "biggest" fragment is conserved
        according to (i) the number of non-Hs atoms and if there is tie then 
        according to (ii) the molecular weight.
        
        :param  mol_in:  RDKit Mol
        :return mol_out: new RDKit Mol having only one connected component
        """
        def count_non_hs_atom(mol):
            ans = 0
            for atm in mol.GetAtoms():
                if atm.GetAtomicNum() != 1:
                    ans += 1
            return ans
        # Remove "other" molecules
        molfrag = GetMolFrags(mol_in, asMols=True, sanitizeFrags=False)
        mol_out = mol_in
        if len(molfrag) > 1:
            accepted_nbr_atm = 0  # flag number of atoms in fragment
            accepted_mw = 0  # flag the molecular weight of the biggest fragment
            for f in molfrag:
                nbr_atm = count_non_hs_atom(f)
                if nbr_atm > accepted_nbr_atm or (nbr_atm == accepted_nbr_atm and MolWt(f) > accepted_mass):
                    accepted_nbr_atm = nbr_atm
                    accepted_mass = MolWt(f)
                    mol_out = f  # keep only the biggest fragment
            cls._copy_properties(mol_in, mol_out)  # save the name and stuff
        return mol_out
        
    @classmethod
    def commute_inchi(cls, mol_in):
        """Convert RDKit compound back and forth to InChi.

        Returns a new compound after the initial one has been converted
        back and forth to InChi.
        
        :param   mol_in:  RDKit Mol
        :return  mol_out: RDKit Mol
        """
        inchi = MolToInchi(mol_in, logLevel=None)  # this is talkative...
        mol_out = MolFromInchi(inchi, sanitize=False, removeHs=False,
                          logLevel=None, treatWarningAsError=False)
        if not mol_out:
            raise ValueError("Failed InChi validity filter.")
        # Copy the properties
        cls._copy_properties(mol_in, mol_out)
        return mol_out

    @classmethod
    def remove_isotope(cls, mol_in):
        """Strip all isotope information.

        Returns a new compound.
        
        :param   mol_in:  RDKit Mol
        :return  mol_out: RDKit Mol
        """
        mol_out = deepcopy(mol_in)  # copy it, just for consistency with other filters
        for atm in mol_out.GetAtoms():
            atm.SetIsotope(0)
        if not mol_out:
            raise ValueError("Failed isotope removing filter.")
        return mol_out

    @staticmethod
    def _rules_rdkit():
        patts = (
            ('[n+;H]', 'n'),  # Imidazoles
            ('[N+;!H0]', 'N'),  # Amines
            ('[$([O-]);!$([O-][#7])]', 'O'),  # Carboxylic acids and alcohols
            ('[S-;X1]', 'S'),  # Thiols
            ('[$([N-;X2]S(=O)=O)]', 'N'),  # Sulfonamides
            ('[$([N-;X2][C,N]=C)]', 'N'),  # Enamines
            ('[n-]', '[nH]'),  # Tetrazoles
            ('[$([S-]=O)]', 'S'),  # Sulfoxides
            ('[$([N-]C=O)]', 'N'),  # Amides
        )
        return [(MolFromSmarts(x), MolFromSmiles(y, False)) for x, y in patts]

    @staticmethod
    def _rules_molvs():
        """Rules to neutralize compounds. Inspired by molvs."""
        ans = {}
        # Neutralizable positive charge (with hydrogens attached)
        # ans["pos_h"] = Chem.MolFromSmarts('[+!H0!$(*~[-])]')
        ans["pos_h"] = MolFromSmarts('[+!H0]')
        # Non-neutralizable positive charge (no hydrogens attached)
        # ans["pos_quat"] = Chem.MolFromSmarts('[+H0!$(*~[-])]')
        # Negative charge, not bonded to a positive charge with no hydrogens
        # ans["neg"] = Chem.MolFromSmarts('[-!$(*~[+H0])]')
        ans["neg"] = MolFromSmarts('[-]')
        # Negative oxygen bonded to [C,P,S]=O, negative aromatic nitrogen?
        # ans["neg_acid"] = Chem.MolFromSmarts('[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]')
        return ans

    @classmethod
    def _neutralise_charge_method1(cls, mol_in, rules=None):
        """Neutralise charges according to a set of predefined rules.

        From:
            http://www.rdkit.org/docs/Cookbook.html#neutralizing-charged-molecules
        """
        # Fallback to default rules if none are provided
        if rules is None:
            fun_rules = cls._rules_rdkit

        # Check if rules are already initialised as an attribute
        if not hasattr(rules, "rules"):
            fun_rules.rules = fun_rules()

        # Apply rules
        # Better to use ReplaceSubstructs than RunReactant: the latter would give
        # several products (or we would need to use HasSubstructMatch anyway).
        for reactant, product in fun_rules.rules:
            while mol_in.HasSubstructMatch(reactant):
                rms = ReplaceSubstructs(mol_in, reactant, product)
                mol_in = rms[0]
        mol_in.UpdatePropertyCache()
        return mol_in

    @classmethod
    def _neutralise_charge_method2(cls, mol_in):
        """Neutralise charges as much as possible playing on hydrogens.

        You should sanitize the compounds after this operation.

        From:
            http://molvs.readthedocs.io/en/latest/_modules/molvs/charge.html
        """
        mol_out = deepcopy(mol_in)  # copy it, just for consistency with other operations
        mol_out.UpdatePropertyCache(strict=False)  # recompute implicit valence
        # Check if rules are already initialised as an attribute
        if not hasattr(cls._rules_molvs, "rules"):
            cls._rules_molvs.rules = cls._rules_molvs()
        # Get atom ids for matches
        p = [x[0] for x in mol_out.GetSubstructMatches(cls._rules_molvs.rules['pos_h'])]
        # q = [x[0] for x in cc.GetSubstructMatches(cls._rules_molvs.rules['pos_quat'])]
        n = [x[0] for x in mol_out.GetSubstructMatches(cls._rules_molvs.rules['neg'])]
        # a = [x[0] for x in cc.GetSubstructMatches(cls._rules_molvs.rules['neg_acid'])]
        # Neutralize negative charges
        # if q:
        #     # Surplus negative charges more than non-neutralizable positive charges
        #     neg_surplus = len(n) - len(q)
        #     if a and neg_surplus > 0:
        #         # zwitterion with more negative charges than quaternary positive centres
        #         while neg_surplus > 0 and a:
        #             # Add hydrogen to first negative acid atom, increase formal charge
        #             # Until quaternary positive == negative total or no more negative acid
        #             atom = cc.GetAtomWithIdx(a.pop(0))
        #             atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        #             atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        #             neg_surplus -= 1
        # Finish of neutralization of negative charges (we don't care for zwitterion)
        for atom in [mol_out.GetAtomWithIdx(x) for x in n]:
            while atom.GetFormalCharge() < 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        # Neutralize positive charges
        for atom in [mol_out.GetAtomWithIdx(x) for x in p]:
            # Remove hydrogen and reduce formal charge until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetTotalNumHs() > 0:
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                if atom.GetNumExplicitHs() > 0:
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
        return mol_out

    @classmethod
    def neutralise_charge(cls, mol_in):
        """Neutralise charges.
        
        :param   mol_in:  RDKit Mol
        :return  mol_out: RDKit Mol
        """
        return cls._neutralise_charge_method1(mol_in)
        # return cls._neutralise_charge_method2(mol_in)

    @classmethod
    def add_hydrogen(cls, mol_in, addCoords=True):
        """Explicit all hydrogens.
        
        :param    mol_in: RDKit Mol
        :param addCoords: Add coordinate to added Hs, bool
        :return  mol_out: RDKit Mol
        """
        return AddHs(mol_in, explicitOnly=False, addCoords=addCoords)

    @classmethod
    def remove_hydrogen(cls, mol_in, addCoords=True):
        """Implicit all hydrogens.
        
        :param    mol_in: RDKit Mol
        :param addCoords: Add coordinate to added Hs, bool
        :return  mol_out: RDKit Mol
        """
        return RemoveHs(mol_in, explicitOnly=False, addCoords=addCoords)

    @classmethod
    def kekulize(cls, mol_in):
        """Kekulize compound.
        
        :param   mol_in:  RDKit Mol
        :return  mol_out: RDKit Mol
        """
        mol_out = deepcopy(mol_in)
        Kekulize(mol_out, clearAromaticFlags=True)
        return mol_out

    @classmethod
    def remove_stereo(cls, mol_in):
        """Wild stereo removal.
        
        Warning: need a back and forth Inchi export/import to normalise tautomer
        
        :param  mol_in:   RDKit mol
        :return mol_out:  RDKit mol
        """
        mol_out = deepcopy(mol_in)
        RemoveStereochemistry(mol_out)
        return mol_out
