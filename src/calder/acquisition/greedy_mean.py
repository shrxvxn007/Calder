"""Highest posterior mean per cost."""

from calder.acquisition.calder import GreedyPolicy


class GreedyMeanPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="greedy_mean", **kwargs)
