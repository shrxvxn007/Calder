"""Upper-confidence bound per cost."""

from calder.acquisition.calder import GreedyPolicy


class UCBPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="ucb", **kwargs)
