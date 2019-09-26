"""
Defines the compound scoring function.
Currently implements toxicity in E. coli, based on data from EcoliTox.
"""

# General use packages
import random
import numpy as np
import sys
import csv
import math
import logging
from rdkit.Chem import DataStructs
from rdkit import Chem

from config import *


class CompoundScoring(object):
    """
    Defines Compound Scorer object.
    """
    logger = logging.getLogger(__name__)
    def __init__(self, scoring_function = None):
        if scoring_function is None:
            pass
        else:
            self.scoring_function = scoring_function

    def __repr__(self):
        """
        Name the used scorer.
        Raises an error is the class is not properly instanciated
        """
        return(self.name)

    def calculate(self, compound):
        score = self.scoring_function(compound)
        return(score)

def pseudo_random(compound):
    score = random.uniform(0, 10)
    return(score)


class ToxicityScoring(CompoundScoring):
    """
    Returns the log toxicity value of a compound.
    The data is stored in a csv file, tab delimited, with columns "name", 'InChI' and "toxicity"
    This can easily be changed to another data with a similar formatting.
    """
    def __init__(self, toxicity_data = "{}/name_structure_toxicity.csv".format(data_path)):
        CompoundScoring.__init__(self)
        self.scoring_function = self.scoring_function()
        self.name = "ToxicityScoring"
        self.fit_model(toxicity_data)

    def calculate_ECFP(self,inchi):
        rdmol = Chem.inchi.MolFromInchi(inchi, sanitize=False)
        # rd_mol = standardize_chemical(rdmol, add_hs=False, heavy = True, rm_stereo=True)
        ECFP= Chem.AllChem.GetMorganFingerprintAsBitVect(rdmol, radius = 2, nBits=1024, useFeatures = False, useChirality = False)
        return(ECFP)

    def select_current_best_model(self, X, y,
                                  models_number = 10,
                                  verbose = False):

        trained_model_list = []
        # Training all models
        for i in range(models_number):
            X_train, y_train = X, y
            other_MLP = MLPRegressor(hidden_layer_sizes  = (10, 100,100, 20), solver ="adam", max_iter=20000,
                                      early_stopping = True, learning_rate = "adaptive")
            other_MLP.fit(X_train, y_train.flatten())
            trained_model_list.append(other_MLP)

            big_MLP = MLPRegressor(hidden_layer_sizes  = (100,100, 20),solver ="adam", max_iter=20000,
                                      early_stopping = True, learning_rate = "adaptive")
            big_MLP.fit(X_train, y_train.flatten())
            trained_model_list.append(big_MLP)


            medium_MLP = MLPRegressor(hidden_layer_sizes  = (40, 10), solver ="adam", max_iter=20000,
                                      early_stopping = True, learning_rate = "adaptive")
            medium_MLP.fit(X_train, y_train.flatten())
            trained_model_list.append(medium_MLP)

            small_MLP = MLPRegressor(hidden_layer_sizes  = (10), solver ="adam", max_iter=20000,
                                      early_stopping = True, learning_rate = "adaptive")
            small_MLP.fit(X_train, y_train.flatten())
            trained_model_list.append(small_MLP)

        # Evaluating all
        all_scores = []
        for i in range(len(trained_model_list)):
            selected_mdoel = trained_model_list[i]
            y_pred = selected_mdoel.predict(X)
            score = sklearn.metrics.r2_score(y, y_pred)
            all_scores.append(score)

        try:
            best_index = all_scores.index(max(all_scores))
            best_score = all_scores[best_index]
        except ValueError:
            best_index = 0
        best_model = trained_model_list[best_index]
        return(best_model, best_score)

    def fit_model(self,toxicity_data):
        y = []
        X = None
        # Loading data
        with open(toxicity_data, "r") as file_hdl:
            reader = csv.DictReader(file_hdl, delimiter = '\t')
            for row in reader:
                y.append(math.log(float(row["toxicity"])))
                arr = np.zeros((1,))
                fp = self.calculate_ECFP(row["InChI"])
                DataStructs.ConvertToNumpyArray(fp, arr)
                arr = np.reshape(arr, (1, 1024))
                if X is None:
                    X = arr
                else:
                    X = np.concatenate((X, arr), axis = 0)
        self.log_loading = "Loaded {} compounds from {}".format(len(y), toxicity_data)
        y = np.array(y)
        # Fitting mdoel:
        best_model, score = self.select_current_best_model(X, y, models_number = 10)
        y_pred = best_model.predict(X)
        score = sklearn.metrics.r2_score(y, y_pred)
        self.log_score = "The toxicity model has a R2 score of {} on itself".format(round(score, 2))
        self.model = best_model

    def scoring_function(self):
        # CODE IT
        def compound_scoring(compound):
            ECFP = compound._get_ECFP()
            arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(ECFP, arr)
            arr = np.reshape(arr, (1, 1024))
            y_pred = self.model.predict(arr)
            return(y_pred)
        return(compound_scoring)


RandomCompoundScorer = CompoundScoring(scoring_function = pseudo_random)
if use_toxicity:
    toxicity_scorer = ToxicityScoring()
