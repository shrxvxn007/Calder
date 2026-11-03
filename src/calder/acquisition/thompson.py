"""Thompson sampling from supplied joint function draws; marginal fallback allowed."""

from calder.acquisition.calder import GreedyPolicy


class ThompsonPolicy(GreedyPolicy):
    def __init__(self, **kwargs: object) -> None:
        super().__init__(name="thompson", **kwargs)
