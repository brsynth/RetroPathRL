"""
Defines the biological scoring function.
Necessitates random for random scoring, and all rule sets for biological scoring.
"""

import random
from rule_sets_examples import *
from rule_sets_similarity import *

class BiologicalScoring(object):
    """
    Defines Biological Scorer object.
    Returns the biological score associated to a reaction rule.
    """
    def __init__(self, scoring_function):
        self.scoring_function = scoring_function
        self.name = "Random"

    def __repr__(self):
        return(self.name)

    def calculate(self, rule):
        score = self.scoring_function(rule)
        return(score)

def pseudo_random(rule):
    score = random.uniform(0, 10)
    return(score)

class BiologicalScoringOrganism(BiologicalScoring):
    """
    Defines Biological Scorer object from an organism with predefined scores.
    Inverted converts a penalty to a score.
    This will be analysed more in depth when biological score will evolve.
    """
    def __init__(self, rules_dictionnary, inverted = False, name = "None"):
        BiologicalScoring.__init__(self, scoring_function = None)
        self.scoring_function = self.assign_from_dict(rules_dictionnary, inverted)
        self.name = name

    def __repr__(self):
        return(self.name)

    def assign_from_dict(self, rules_dictionnary, inverted):
        rules_dictionnary = rules_dictionnary
        def simple_assign_inside(rule):
            score = rules_dictionnary[rule]["biological_score"]
            # Inverted if to use penalties instead of scors.
            # if inverted:
            #     try:
            #         return(1/score)
            #     except ZeroDivisionError:
            #         return(33)
            # else:
            return(score)
        return(simple_assign_inside)


RandomBiologicalScorer = BiologicalScoring(scoring_function = pseudo_random)
BiologicalFullScoringRetroH = BiologicalScoringOrganism(rules_dictionnary= full_rules_retro_H, name = "full_rules_retro_H")
BiologicalFullScoringFwdH = BiologicalScoringOrganism(rules_dictionnary= full_rules_forward_H, name = "full_rules_forward_H")
BiologicalFullScoringRetroNoH = BiologicalScoringOrganism(rules_dictionnary= full_rules_retro_no_H, name = "full_rules_retro_no_H")
BiologicalFullScoringFwdNoH = BiologicalScoringOrganism(rules_dictionnary= full_rules_forward_no_H, name = "full_rules_forward_no_H")

full_H = full_rules_retro_H
full_H.update(full_rules_forward_H)
BiologicalFullScoringH = BiologicalScoringOrganism(rules_dictionnary= full_H, name = "full_rules_retro_H")
