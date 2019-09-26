"""
Defines the UCT (Upper Confidence Tree) policies.
It it the formula that allows for balancing between exploration and exploitation when selecting children in the Tree.
Implements a number of different policies.
Policies are Subclasses of UCT_policy Class.
They need to have a attribute function that does the calculation. See examples if you want to develop your own.
"""

from math import sqrt, log


class UCT_policy(object):
    """
    Defines UCT_policies objects.
    They take a node and return the best child according to this policy.
    Only subclasses of this object can work as there is no default calculation function.
    """
    def __init__(self, parameters = {"UCTK": 2}, policy_type = 'Classical', function = None):
        self.parameters = parameters
        self.policy_type = policy_type

    def calculate(self, node, top_n = 1):
        s = sorted(node.children, key = lambda c: self.function(c, parent_visits = node.visits))
        s = s[-top_n]
        return s

    def __str__(self):
        return("Policy type: {} \nFormula: {}".format(self.policy_type, self.formula))


class Classical_UCT(UCT_policy):
    """
    This class implements the most basic UCT functions.
    Only uses number of visits as a criteria.
    It is the Classical UCT formula where no additionnal expert knowledge is inputed.
    """
    def __init__(self, parameters = {"UCTK": 1000}):
        UCT_policy.__init__(self, policy_type = "Classical")
        self.parameters = parameters
        self.formula = "mean_score + sqrt({}*log(N + 1)/(n+1))".format(parameters["UCTK"])
        self.function = self.simple_UCT_formula(self.parameters)

    def simple_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        def simple_formula_inside(c, parent_visits):
            value = c.average_score + sqrt(UCTK*log(parent_visits +1)/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Classical_UCT_RAVE(UCT_policy):
    """
    This class implements UCT based on visit count and RAVE.
    RAVE stands for Rapid Action Value Estimation:
    - it adds another score based on usage of identical moves elsewhere in the Tree
    - this is ponderated by the number of visits: as visits increase, the actual score of the node becomes more important than this initila estimation.
    """
    def __init__(self, parameters = {"UCTK": 1000, "k_rave": 100}):
        UCT_policy.__init__(self, policy_type = "Classical_RAVE")
        self.parameters = parameters
        self.formula = "(1-b) mean_score + b rave_score + sqrt({}*log(N + 1)/(n+1)) with b = sqrt({}/(3N + {}))".format(parameters["UCTK"], parameters["k_rave"], parameters["k_rave"])
        self.function = self.RAVE_formula(parameters = self.parameters)

    def RAVE_formula(self, parameters):
        UCTK = parameters["UCTK"]
        k_rave = parameters["k_rave"]
        def simple_formula_inside(c, parent_visits):
            b = sqrt(k_rave/(3*parent_visits + k_rave))
            value = c.average_score *(1-b) + b * c.move.RAVE_average_score + sqrt(UCTK*log(parent_visits +1)/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Classical_UCT_with_bias(UCT_policy):
    """
    This class implements UCT based on visits and progressive bias.
    Progressive bias works by
    - giving an initial value to a node (based on expert knowledge for example)
    - this importance decreases as the node gets visited and this initial estimation's importance decreases in favor of actual rollouts.
    """
    def __init__(self, parameters = {"UCTK": 1000, "bias_k": 1}):
        UCT_policy.__init__(self, policy_type = "Classical")
        self.parameters = parameters
        self.formula = "mean_score + sqrt({}*log(N + 1)/(n+1)) + {} *  progressive_bias/(n+1)".format(parameters["UCTK"], parameters["bias_k"])
        self.function = self.simple_UCT_formula(self.parameters)

    def simple_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        bias_k = parameters["bias_k"]
        def simple_formula_inside(c, parent_visits):
            value = c.average_score + sqrt(UCTK*log(parent_visits +1)/(c.visits + 1))+ bias_k * c.progressive_bias/(c.visits + 1)
            return(value)
        return(simple_formula_inside)

class Nature_UCT(UCT_policy):
    """
    This class implements the formula used in the following Nature paper :(https://doi.org/10.1038/nature25978)
    Planning chemical syntheses with deep neural networks and symbolic AI
    It is identical to the Chemical Scoring UCT (Chemical_UCT_1)
    """
    def __init__(self, parameters = {"UCTK": 3}):
        UCT_policy.__init__(self, policy_type = "Nature Symbolic IA")
        self.parameters = parameters
        self.formula = "mean_score + {} * P * sqrt(N/(n+1))".format(parameters["UCTK"])
        self.function = self.Nature_UCT_formula(self.parameters)

    def Nature_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        def simple_formula_inside(c, parent_visits):
            chem_P = c.move.chemical_score
            value = c.average_score + UCTK * chem_P *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Biochemical_UCT_1(UCT_policy):
    """
    This class implements a simple biochemical score UCT.
    The selection is guided by a product of chemical and biological score.
    """
    def __init__(self, parameters = {"UCTK": 3}):
        UCT_policy.__init__(self, policy_type = "Biochemical multiplication")
        self.parameters = parameters
        self.formula = "mean_score + {} * P_c * B * sqrt(N/(n+1))".format(parameters["UCTK"])
        self.function = self.Biochemical_UCT_formula(self.parameters)

    def Biochemical_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        def simple_formula_inside(c, parent_visits):
            chem_P = c.move.chemical_score
            b_score = c.move.biological_score
            value = c.average_score + UCTK * chem_P * b_score *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Biological_UCT_1(UCT_policy):
    """
    This class implements a simple biological score UCT.
    The selection is guided by Biological score only.
    """
    def __init__(self, parameters = {"UCTK": 3}):
        UCT_policy.__init__(self, policy_type = "Biological score only")
        self.parameters = parameters
        self.formula = "mean_score + {} * B * sqrt(N/(n+1))".format(parameters["UCTK"])
        self.function = self.Biological_UCT_formula(self.parameters)

    def Biological_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        def simple_formula_inside(c, parent_visits):
            b_score = c.move.biological_score
            value = c.average_score + UCTK * b_score *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Chemical_UCT_1(UCT_policy):
    """
    This class implements a simple chemical score UCT.
    The selection is guided by Chemical score only.
    """
    def __init__(self, parameters = {"UCTK": 3}):
        UCT_policy.__init__(self, policy_type = "Chemical multiplication")
        self.parameters = parameters
        self.formula = "mean_score + {} * P_c * sqrt(N/(n+1))".format(parameters["UCTK"])
        self.function = self.Chemical_UCT_formula(self.parameters)

    # @staticmethod
    def Chemical_UCT_formula(self, parameters):
        UCTK = parameters["UCTK"]
        def simple_formula_inside(c, parent_visits):
            chem_P = c.move.chemical_score
            value = c.average_score + UCTK * chem_P *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Biochemical_UCT_1_with_RAVE(UCT_policy):
    """
    This class implements a biochemical score UCT with RAVE augmentation.
    RAVE stands for Rapid Action Value Estimation:
    - it adds another score based on usage of identical moves elsewhere in the Tree
    - this is ponderated by the number of visits: as visits increase, the actual score of the node becomes more important than this initila estimation.
    """
    def __init__(self, parameters = {"UCTK": 3, "k_rave": 100}):
        UCT_policy.__init__(self, policy_type = "Biochemical multiplication with RAVE")
        self.parameters = parameters
        self.formula = "(1-b) mean_score + b rave_score + {} * P_c * B * sqrt(N/(n+1)) with b = sqrt({}/(3N + {})".format(parameters["UCTK"], parameters["k_rave"], parameters["k_rave"])
        self.function = self.Biochemical_UCT_RAVE_formula(self.parameters)

    def Biochemical_UCT_RAVE_formula(self, parameters):
        UCTK = parameters["UCTK"]
        k_rave = parameters["k_rave"]
        def simple_formula_inside(c, parent_visits):
            b = sqrt(k_rave/(3*parent_visits + k_rave))
            b_score = c.move.biological_score
            chem_P = c.move.chemical_score
            value = c.average_score * (1-b) + b * c.move.RAVE_average_score + UCTK * chem_P * b_score *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Biochemical_UCT_with_progressive_bias(UCT_policy):
    """
    This class implements a biochemical score UCT and progressive bias.
    Progressive bias works by
    - giving an initial value to a node (based on expert knowledge for example)
    - this importance decreases as the node gets visited and this initial estimation's importance decreases in favor of actual rollouts.
    """
    def __init__(self, parameters = {"UCTK": 3, "bias_k": 1}):
        UCT_policy.__init__(self, policy_type = "Biochemical with progressive bias")
        self.parameters = parameters
        self.formula = "mean_score + {} * bias/(n+1) + {} * P_c * B * sqrt(N/(n+1))".format(parameters["bias_k"], parameters["UCTK"])
        self.function = self.Biochemical_UCT_with_bias_formula(parameters)

    # @staticmethod
    def Biochemical_UCT_with_bias_formula(self, parameters):
        UCTK = parameters["UCTK"]
        bias_k = parameters["bias_k"]
        def simple_formula_inside(c, parent_visits):
            chem_P = c.move.chemical_score
            b_score = c.move.biological_score
            bias = c.progressive_bias
            value = c.average_score + bias_k * bias/(c.visits +1)  + UCTK * chem_P * b_score *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)

class Biochemical_UCT_with_toxicity(UCT_policy):
    """
    This class implements a biochemical score UCT combined with toxicity bias.
    the formula is identical to the Biochemical_UCT_with_progressive_bias, the bias being the node's toxicity.
    """
    def __init__(self, parameters = {"UCTK": 3, "bias_k": 1}):
        UCT_policy.__init__(self, policy_type = "Biochemical with toxicity")
        self.parameters = parameters
        self.formula = "mean_score + {} * toxicity/(n+1) + {} * P_c * B * sqrt(N/(n+1))".format(parameters["bias_k"], parameters["UCTK"])
        self.function = self.Biochemical_UCT_with_toxicity_formula(parameters)

    def Biochemical_UCT_with_toxicity_formula(self, parameters):
        UCTK = parameters["UCTK"]
        bias_k = parameters["bias_k"]
        def simple_formula_inside(c, parent_visits):
            chem_P = c.move.chemical_score
            b_score = c.move.biological_score
            toxicity = c.toxicity
            value = c.average_score + bias_k * toxicity/(c.visits +1)  + UCTK * chem_P * b_score *sqrt(parent_visits/(c.visits + 1))
            return(value)
        return(simple_formula_inside)
