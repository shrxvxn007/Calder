"""Optional read-only view of saved retrospective outputs."""

import sys
from pathlib import Path

import pandas as pd
import streamlit as st

root = Path(sys.argv[1])
st.title("CALDER saved experiment explorer")
st.caption(
    "Retrospective simulated selections. Observed outcomes, predictions and epistemic uncertainty are distinct. No retraining occurs."
)
rounds = pd.read_csv(root / "rounds.csv")
selected = pd.read_csv(root / "selections.csv")
scenario = st.selectbox("Dataset / scenario", sorted(rounds.scenario.unique()))
policy = st.selectbox("Policy", sorted(rounds.policy.unique()))
seed = st.selectbox("Seed", sorted(rounds.seed.unique()))
q = st.selectbox("Pilot-derived threshold quantile", sorted(rounds.threshold_quantile.unique()))
mask = (
    (rounds.scenario == scenario)
    & (rounds.policy == policy)
    & (rounds.seed == seed)
    & (rounds.threshold_quantile == q)
)
st.line_chart(rounds[mask].set_index("additional_budget")[["confirmed_hits", "group_coverage"]])
choices = rounds[mask]["round"].tolist()
current = st.select_slider("Round", options=choices)
mask = (
    (selected.scenario == scenario)
    & (selected.policy == policy)
    & (selected.seed == seed)
    & (selected.threshold_quantile == q)
    & (selected["round"] == current)
)
st.dataframe(selected[mask])
