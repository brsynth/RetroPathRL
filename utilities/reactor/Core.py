"""
Core code for firing rules
"""


class RuleMatchError(Exception):
    """Raised when something went wrong when matching a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-MATCH-ERROR: {}".format(self._msg)


class RuleFireError(Exception):
    """Raised when something went wrong when firing a rule."""

    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return "RULE-FIRE-ERROR: {}".format(self._msg)


class RuleBurnerCore(object):
    """Apply one rule on one chemical.""" 

    def __init__(self, rd_rule, rd_mol):
        """Apply one rule on one chemical.
        
        Notice: no standardization is made on inputed chemicals and rules.
        
        :param  rd_rule:    RDKit reaction object, reactio rule to apply
        :param  rd_mol:     RDKit mol object, chemical
        :param  timeout:    str, Reaction rule SMARTS
        """
        # Internal settings
        USE_CHIRALITY_IN_MATCH = False  # Default value anyway substrucre matching
        # Input
        self._rd_rule = rd_rule
        self._rd_mol = rd_mol
    
    def match(self):
        """Check if left reaction side match the chemical.
        
        returns:    bool, True if there is a match, else False
        """
        try:
            for reactant in self._rd_rule.GetReactants():
                if self._rd_mol.HasSubstructMatch(reactant, ):
                    return True
            return False
        except Exception as e:
            raise RuleMatchError(e) from e
        
    def fire(self):
        """Fire the rule on the chemical.
        
        returns:    tuple of tuple, list of results for each possible application.
        """
        try:
            return self._rd_rule.RunReactants((self._rd_mol,))
        except Exception as e:
            raise RuleFireError(e) from e
