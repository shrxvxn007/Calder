"""Independently implemented Probability-of-Hit baseline."""

from calder.acquisition.calder import GreedyPolicy


class ProbabilityOfHitPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="probability_of_hit", **kwargs)
