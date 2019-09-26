
"""
Defines the Rollout policies.
Usage is : move = RolloutPolicy.select_best_move(available_moves)
Remarks:
- various policies have been tested on toy examples on a Jupyter notebook during implementation
"""

from math import sqrt, log
import random

class Rollout_policy(object):
    """
    Defines rollout policy.
    From a list of moves, select the one that should be used for rollout.
    This is the base object, subclasses necessitate a policy function.
    """
    def __init__(self, policy_type, description = "Default Rollout Policy"):
        self.policy_type = policy_type
        self.description = description

    def select_best_move(self, available_moves):
        try:
            move = self.policy(available_moves)
            return(move)
        except IndexError:
            return(None)

    def __str__(self):
        return("Policy type: {} \nDescription: {}".format(self.policy_type, self.description))

class Rollout_policy_first(Rollout_policy):
    """
    Defines rollout policy.
    Always returns the first element: first compound, first rule
    """
    def __init__(self):
        description = "Always select the first compound_rule combination"
        Rollout_policy.__init__(self, policy_type = "First found combination", description = description)
        self.name = "Rollout_policy_first"
        self.policy = self.policy()

    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            move = available_moves[0]
            return(move)
        return(select_best_inside)

class Rollout_policy_chemical_best(Rollout_policy):
    """
    Defines rollout policy.
    Always returns the best chemical move
    """
    def __init__(self):
        description = "Always select the move with the highest chemical score"
        Rollout_policy.__init__(self, policy_type = "Best Chemical", description = description)
        self.policy = self.best_chemical_policy()
        self.name = "Rollout_policy_chemical_best"

    def best_chemical_policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            current_best = available_moves[0]
            current_best_score = current_best.chemical_score
            for element in available_moves:
                chemical_score = element.chemical_score
                if chemical_score > current_best_score:
                    current_best_score = chemical_score
                    current_best = element
            return(current_best)
        return(select_best_inside)

class Rollout_policy_biological_best(Rollout_policy):
    """
    Defines rollout policy.
    Always returns the best biological move
    """
    def __init__(self):
        description = "Always select the move with the highest biological score"
        Rollout_policy.__init__(self, policy_type = "Best Biological", description = description)
        self.policy = self.best_biological_policy()
        self.name = "Rollout_policy_biological_best"

    def best_biological_policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            current_best = available_moves[0]
            current_best_score = current_best.biological_score
            for element in available_moves:
                biological_score = current_best_score = element.biological_score
                if biological_score > current_best_score:
                    current_best_score = biological_score
                    current_best = element
            return(current_best)
        return(select_best_inside)

class Rollout_policy_biochemical_addition_best(Rollout_policy):
    """
    Defines rollout policy.
    Always returns the best biochemical (addition of scores) move
    """
    def __init__(self):
        description = "Select the highest Biochemical addition score"
        Rollout_policy.__init__(self, policy_type = "Best Biochemical addition", description = description)
        self.policy = self.best_biochemical_policy()
        self.name = "Rollout_policy_biochemical_addition_best"

    def best_biochemical_policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            current_best = available_moves[0]
            current_best_score = current_best.biological_score + current_best.chemical_score
            for element in available_moves:
                biological_score = element.biological_score
                chemical_score = element.chemical_score
                if biological_score + chemical_score > current_best_score:
                    current_best_score = biological_score + chemical_score
                    current_best = element
            return(current_best)
        return(select_best_inside)

class Rollout_policy_biochemical_multiplication_best(Rollout_policy):
    """
    Defines rollout policy.
    Always returns the best biochemical (multiplication of scores) move
    """
    def __init__(self):
        description = "Select the highest Biochemical multiplication score"
        Rollout_policy.__init__(self, policy_type = "Best Biochemical multiplication", description = description)
        self.policy = self.best_biochemical_policy()
        self.name = "Rollout_policy_biochemical_multiplication_best"

    def best_biochemical_policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            current_best = available_moves[0]
            current_best_score = current_best.biological_score * current_best.chemical_score
            for element in available_moves:
                biological_score = element.biological_score
                chemical_score = element.chemical_score
                if biological_score * chemical_score > current_best_score:
                    current_best_score = biological_score * chemical_score
                    current_best = element
            return(current_best)
        return(select_best_inside)

class Rollout_policy_random_uniform(Rollout_policy):
    """
    Random sampling of the move amongst available moves
    """
    def __init__(self):
        description = "Random selection - no scoring involved"
        Rollout_policy.__init__(self, policy_type = "Random sampling", description = description)
        self.policy = self.policy()
        self.name = "Rollout_policy_random_uniform"

    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            index = random.randrange(0, len(available_moves))
            move = available_moves[index]
            return(move)
        return(select_best_inside)

class Rollout_policy_random_uniform_on_chem_score(Rollout_policy):
    """
    Random sampling of the move amongst available moves, weighted by chemical score
    """
    def __init__(self):
        description = "Random selection - uniform sampling from chemical weights"
        Rollout_policy.__init__(self, policy_type = "Chemical uniform sampling", description = description)
        self.policy = self.policy()
        self.name = "Rollout_policy_random_uniform_on_chem_score"

    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            pop, cum, cum_w = [], [], 0
            for move in available_moves:
                pop.append(move)
                cum_w = cum_w + move.chemical_score
                cum.append(cum_w)
            move = random.choices(pop, cum_weights=cum, k=1)[0]
            return(move)
        return(select_best_inside)

class Rollout_policy_random_uniform_on_bio_score(Rollout_policy):
    """
    Random sampling of the move amongst available moves, weighted by biological score
    """
    def __init__(self):
        description = "Random selection - uniform sampling from biological weights"
        Rollout_policy.__init__(self, policy_type = "Biological uniform sampling", description = description)
        self.policy = self.policy()
        self.name = "Rollout_policy_random_uniform_on_bio_score"
    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            pop, cum, cum_w = [], [], 0

            for move in available_moves:
                pop.append(move)
                cum_w = cum_w + move.biological_score
                cum.append(cum_w)
            move = random.choices(pop, cum_weights=cum, k=1)[0]
            return(move)
        return(select_best_inside)

class Rollout_policy_random_uniform_on_biochemical_addition_score(Rollout_policy):
    """
    Random sampling of the move amongst available moves, weighted by biochemical (addition) score
    """
    def __init__(self):
        description = "Random selection - uniform sampling from added biochemical weights"
        Rollout_policy.__init__(self, policy_type = "Biochemical addition uniform sampling", description = description)
        self.policy = self.policy()
        self.name = "Rollout_policy_random_uniform_on_biochemical_addition_score"

    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            pop, cum, cum_w = [], [], 0

            for move in available_moves:
                pop.append(move)
                cum_w = cum_w + move.biological_score + move.chemical_score
                cum.append(cum_w)
            move = random.choices(pop, cum_weights=cum, k=1)[0]
            return(move)
        return(select_best_inside)

class Rollout_policy_random_uniform_on_biochemical_multiplication_score(Rollout_policy):
    """
    Random sampling of the move amongst available moves, weighted by biochemical (multiplication) score
    """
    def __init__(self):
        description = "Random selection - uniform sampling from multiplied biochemical weights"
        Rollout_policy.__init__(self, policy_type = "Biochemical uniform sampling", description = description)
        self.policy = self.policy()
        self.name = "Rollout_policy_random_uniform_on_biochemical_multiplication_score"

    def policy(self):
        # CODE IT
        def select_best_inside(available_moves):
            pop, cum, cum_w = [], [], 0

            for move in available_moves:
                pop.append(move)
                cum_w = cum_w + move.biological_score * move.chemical_score
                cum.append(cum_w)
            move = random.choices(pop, cum_weights=cum, k=1)[0]
            return(move)
        return(select_best_inside)
