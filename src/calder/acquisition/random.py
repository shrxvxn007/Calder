"""Uniform random ordering with budget feasibility filtering."""

from calder.acquisition.calder import GreedyPolicy


class RandomPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="random", **kwargs)
