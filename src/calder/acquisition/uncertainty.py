"""Largest epistemic standard deviation per cost."""

from calder.acquisition.calder import GreedyPolicy


class UncertaintyPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="uncertainty", **kwargs)
