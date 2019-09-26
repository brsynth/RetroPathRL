"""
Defines the pathway scoring functions.
Can take as inputs both Pathway objects and json dictionnaries exported from Pathways.
"""

import random
import numpy as np
import json
import os
# RP3 - specific objects
from pathway import Pathway


def geo_mean(iterable):
    a = np.array(iterable)
    return a.prod()**(1.0/len(a))

# def geo_mean_overflow(iterable):
#     a = np.log(iterable)
#     return np.exp(a.sum()/len(a))

class PathwayScoring(object):
    """
    Defines Pathway Scorer object.
    """
    def __init__(self, scoring_function = None, scoring_json_function = None):
        if scoring_function is None:
            pass
        else:
            self.scoring_function = scoring_function
        if scoring_json_function is None:
            pass
        else:
            self.scoring_json_function = scoring_json_function

    def __repr__(self):
        """
        Name the used scorer.
        Raises an error is the class is not properly instantiated
        """
        return(self.name)

    def calculate(self, pathway):
        score = self.scoring_function(pathway)
        return(score)

    def calculate_json(self, pathway):
        score = self.scoring_json_function(pathway)
        return(score)

def pseudo_random(pathway):
    score = random.uniform(0, 10)
    return(score)

class ConstantPathwayScoring(PathwayScoring):
    """
    Returns a constant reward, whichever the pathway.
    """
    def __init__(self, reward = 10):
        PathwayScoring.__init__(self)
        self.reward = reward
        self.scoring_function = self.scoring_function()
        self.scoring_json_function = self.scoring_json_function()
        self.name = "ConstantPathwayScoring of {}".format(reward)

    def set_reward(self,reward):
        # For changing the reward of the object
        self.reward = reward
        self.scoring_function = self.scoring_function()
        self.scoring_json_function = self.scoring_json_function()

    def scoring_function(self):
        def pathway_scoring(pathway):
            return(self.reward)
        return(pathway_scoring)

    def scoring_json_function(self):
        def pathway_scoring(pathway):
            return(self.reward)
        return(pathway_scoring)

class BiologicalPathwayScoring(PathwayScoring):
    """
    Returns the geometric mean of biological scores in the Pathway.
    """
    def __init__(self):
        PathwayScoring.__init__(self)
        self.scoring_function = self.scoring_function()
        self.scoring_json_function = self.scoring_json_function()
        self.name = "BiologicalPathwayScoring"

    def scoring_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway.nodes_transformations:
                scores.append(move["data"]["Score"])
            return(geo_mean(scores))
        return(pathway_scoring)

    def scoring_json_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway["elements"]["nodes"]:
                if move["data"]["type"] == "reaction":
                    scores.append(move["data"]["Score"])
            return(geo_mean(scores))
        return(pathway_scoring)

class ChemicalPathwayScoring(PathwayScoring):
    """
    Returns the geometric mean of chemical scores in the Pathway.
    """
    def __init__(self):
        PathwayScoring.__init__(self)
        self.scoring_function = self.scoring_function()
        self.scoring_json_function = self.scoring_json_function()
        self.name = "ChemicalPathwayScoring"

    def scoring_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway.nodes_transformations:
                scores.append(move["data"]["ChemicalScore"])
            return(geo_mean(scores))
        return(pathway_scoring)

    def scoring_json_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway["elements"]["nodes"]:
                if move["data"]["type"] == "reaction":
                    scores.append(move["data"]["ChemicalScore"])
            return(geo_mean(scores))
        return(pathway_scoring)

class BiochemicalPathwayScoring(PathwayScoring):
    """
    Returns the geometric mean of biochemical scores in the Pathway.
    """
    def __init__(self):
        PathwayScoring.__init__(self)
        self.scoring_function = self.scoring_function()
        self.scoring_json_function = self.scoring_json_function()
        self.name = "ChemicalPathwayScoring"

    def scoring_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway.nodes_transformations:
                scores.append(move["data"]["ChemicalScore"] * move["data"]["Score"])
            return(geo_mean(scores))
        return(pathway_scoring)

    def scoring_json_function(self):
        def pathway_scoring(pathway):
            scores = []
            for move in pathway["elements"]["nodes"]:
                if move["data"]["type"] == "reaction":
                    scores.append(move["data"]["Score"] * move["data"]["ChemicalScore"])
            return(geo_mean(scores))
        return(pathway_scoring)

RandomPathwayScorer = PathwayScoring(scoring_function = pseudo_random)
constant_pathway_scoring = ConstantPathwayScoring(reward = 10)
null_pathway_scoring = ConstantPathwayScoring(reward = 0)
biological_pathway_scoring = BiologicalPathwayScoring()
chemical_pathway_scoring = ChemicalPathwayScoring()
biochemical_pathway_scoring = BiochemicalPathwayScoring()
