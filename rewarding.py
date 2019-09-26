"""
Defines the possible rewards for rollout.
Can be augmented for more complex policies using simialr scheme as Rollout or UCT policies.
Is defined through CLI in the Tree script.
"""

class RolloutRewards(object):
    """
    Defines penalty and rewards for the rollout if it's in the chasis.
    """
    def __init__(self, penalty, full_state_reward):
        self.penalty = penalty
        self.full_state_reward = full_state_reward

    def __repr__(self):
        """Reward representation is its values"""
        return("Penalty is {} and full state reward is {}".format(self.penalty, self.full_state_reward))

Basic_Rollout_Reward = RolloutRewards(penalty = -1, full_state_reward = 2)
