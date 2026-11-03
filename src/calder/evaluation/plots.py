"""Figures generated exclusively from saved experiment tables."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from calder.evaluation.calibration import reliability_bins
from calder.evaluation.statistics import bootstrap_interval

STRATA = ("dataset", "scenario", "threshold_quantile")


def _slug(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "all"


def _strata(frame: pd.DataFrame):
    keys = [k for k in STRATA if k in frame]
    if not keys:
        yield "all", frame
        return
    for group, part in frame.groupby(keys, dropna=False, sort=True):
        values = group if isinstance(group, tuple) else (group,)
        yield "__".join(f"{k}-{_slug(v)}" for k, v in zip(keys, values, strict=True)), part


def _save(fig, folder: Path, name: str) -> str:
    fig.tight_layout()
    path = folder / f"{name}.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return str(path.name)


def _curve_plot(part: pd.DataFrame, metric: str, ylabel: str, title: str):
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    budget = "additional_budget" if "additional_budget" in part else "budget_spent"
    for policy, subset in part.groupby("policy", sort=True):
        seed_curves = []
        # Common support only; never extrapolate shorter policy runs.
        starts, stops = [], []
        for _, run in subset.groupby("seed", sort=True):
            run = run.sort_values(budget).drop_duplicates(budget, keep="last")
            starts.append(run[budget].min())
            stops.append(run[budget].max())
        if not starts:
            continue
        grid = np.unique(subset[budget].to_numpy(float))
        grid = grid[(grid >= max(starts)) & (grid <= min(stops))]
        for _, run in subset.groupby("seed", sort=True):
            run = run.sort_values(budget).drop_duplicates(budget, keep="last")
            idx = np.searchsorted(run[budget], grid, side="right") - 1
            seed_curves.append(run[metric].to_numpy(float)[idx])
        values = np.asarray(seed_curves)
        if not grid.size:
            continue
        means, lows, highs = [], [], []
        for column in values.T:
            mean, lo, hi = bootstrap_interval(column, n_resamples=1000)
            means.append(mean)
            lows.append(lo)
            highs.append(hi)
        line = ax.step(grid, means, where="post", label=str(policy))[0]
        ax.fill_between(grid, lows, highs, step="post", color=line.get_color(), alpha=0.12)
    ax.set(xlabel="Additional experimental cost beyond pilot", ylabel=ylabel, title=title)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    ax.grid(alpha=0.18)
    return fig


def plot_saved_results(run_dir: Path | str) -> dict[str, dict[str, object]]:
    """Produce available figures and an explicit status for every requested figure."""
    root = Path(run_dir)
    folder = root / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    status = {
        f"figure_{k}": {
            "status": "unavailable",
            "reason": "Required saved inputs absent",
            "files": [],
        }
        for k in range(1, 9)
    }
    rounds_path = root / "rounds.csv"
    rounds = pd.read_csv(rounds_path) if rounds_path.exists() else pd.DataFrame()
    preq_path = root / "prequential.csv"
    preq = (
        pd.read_csv(preq_path)
        if preq_path.exists() and preq_path.stat().st_size > 1
        else pd.DataFrame()
    )
    endpoints_path = root / "seed_level_results.csv"
    ends = pd.read_csv(endpoints_path) if endpoints_path.exists() else pd.DataFrame()
    with plt.rc_context(
        {"font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False}
    ):
        if not rounds.empty:
            for name, part in _strata(rounds):
                for number, metric, label in [
                    (1, "confirmed_hits", "Cumulative confirmed measured hits"),
                    (2, "group_coverage", "Distinct fixed groups with confirmed hits"),
                ]:
                    if metric in part:
                        fig = _curve_plot(part, metric, label, name.replace("__", " | "))
                        status[f"figure_{number}"]["files"].append(
                            _save(fig, folder, f"figure_{number}_{name}")
                        )
        if not preq.empty and {"probability", "is_hit", "raw_probability"}.issubset(preq.columns):
            for name, part in _strata(preq):
                fig, ax = plt.subplots(figsize=(6.5, 5))
                # The same adaptively selected examples supply both curves.
                # Prefer CALDER examples to avoid conflating policy populations.
                selected = (
                    part[part.policy == "calder"]
                    if "policy" in part and (part.policy == "calder").any()
                    else part
                )
                target = (
                    "latent_label"
                    if "latent_label" in selected and selected.latent_label.notna().all()
                    else "is_hit"
                )
                calibrated_column = (
                    "calibrated_probability"
                    if "calibrated_probability" in selected
                    else "probability"
                )
                for column, label in [
                    ("raw_probability", "Raw latent hit probability"),
                    (calibrated_column, "Pilot measured-label calibration"),
                ]:
                    table = reliability_bins(selected[column], selected[target].astype(int))
                    used = table[table["count"] > 0]
                    ax.plot(used.mean_probability, used.hit_frequency, "o-", label=label)
                ax.plot([0, 1], [0, 1], ":", color="0.5")
                ax.set(
                    xlim=(0, 1),
                    ylim=(0, 1),
                    xlabel="Mean predicted hit probability",
                    ylabel="Empirical hit frequency",
                    title=f"{name}\nPrequential, selected examples; target={target}",
                )
                ax.legend(fontsize=8, frameon=False)
                status["figure_3"]["files"].append(_save(fig, folder, f"figure_3_{name}"))
        if not ends.empty and {"confirmed_hits", "group_coverage", "policy"}.issubset(ends.columns):
            for name, part in _strata(ends):
                fig, ax = plt.subplots(figsize=(6.6, 4.8))
                for policy, run in part.groupby("policy", sort=True):
                    x, xlo, xhi = bootstrap_interval(run.confirmed_hits)
                    y, ylo, yhi = bootstrap_interval(run.group_coverage)
                    ax.errorbar(
                        x,
                        y,
                        xerr=[[x - xlo], [xhi - x]],
                        yerr=[[y - ylo], [yhi - y]],
                        fmt="o",
                        capsize=2,
                        label=str(policy),
                    )
                ax.set(
                    xlabel="Final confirmed measured hits",
                    ylabel="Final distinct groups with confirmed hits",
                    title=name.replace("__", " | "),
                )
                ax.legend(fontsize=7, ncol=2, frameon=False)
                status["figure_4"]["files"].append(_save(fig, folder, f"figure_4_{name}"))
        ablation_path = root / "ablations.csv"
        ablations = (
            pd.read_csv(ablation_path)
            if ablation_path.exists() and ablation_path.stat().st_size > 1
            else pd.DataFrame()
        )
        if not ablations.empty and "confirmed_hits" in ablations:
            variant = next((k for k in ("ablation", "variant", "policy") if k in ablations), None)
            if variant:
                for name, part in _strata(ablations):
                    means = part.groupby(variant).confirmed_hits.mean().sort_index()
                    fig, ax = plt.subplots(figsize=(8, 4.8))
                    ax.barh(means.index.astype(str), means.values)
                    ax.set(
                        xlabel="Mean final confirmed measured hits",
                        title=f"Controlled ablations | {name}",
                    )
                    status["figure_5"]["files"].append(_save(fig, folder, f"figure_5_{name}"))
        if not ends.empty and "dataset" in ends:
            bio = ends[ends.dataset.astype(str).str.contains("replogle22", regex=False)]
            if bio.dataset.nunique() >= 2:
                context_keys = [k for k in ("threshold_quantile", "model") if k in bio]
                context_parts = (
                    bio.groupby(context_keys, dropna=False, sort=True)
                    if context_keys
                    else [((), bio)]
                )
                for context, matching in context_parts:
                    if matching.dataset.nunique() < 2:
                        continue
                    name = _slug(context)
                    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
                    for ax, metric in zip(axes, ("confirmed_hits", "group_coverage"), strict=True):
                        pivot = (
                            matching.groupby(["policy", "dataset"])[metric]
                            .mean()
                            .unstack("dataset")
                        )
                        pivot.plot.bar(ax=ax)
                        ax.set_ylabel(metric.replace("_", " "))
                        ax.tick_params(axis="x", labelrotation=65, labelsize=7)
                        ax.set_xlabel("")
                    fig.suptitle(
                        f"Independent context reruns, frozen method hyperparameters | {name}"
                    )
                    status["figure_6"]["files"].append(_save(fig, folder, f"figure_6_{name}"))
            else:
                status["figure_6"]["reason"] = (
                    "Both K562 and RPE1 measured benchmark endpoints are required"
                )
        if not preq.empty and {"mean", "epistemic_variance", "observed_effect"}.issubset(
            preq.columns
        ):
            for name, part in _strata(preq):
                target = next(
                    (
                        k
                        for k in ("latent_effect", "true_log_effect")
                        if k in part and part[k].notna().all()
                    ),
                    "observed_effect",
                )
                fig, ax = plt.subplots(figsize=(6.6, 4.7))
                for policy, rows in part.groupby("policy", sort=True):
                    ax.scatter(
                        np.sqrt(np.maximum(rows.epistemic_variance, 0)),
                        np.abs(rows["mean"] - rows[target]),
                        s=9,
                        alpha=0.35,
                        label=str(policy),
                    )
                ax.set(
                    xlabel="Pre-query epistemic standard deviation (log effect)",
                    ylabel=f"Absolute prediction error vs {target}",
                    title=name.replace("__", " | "),
                )
                ax.legend(fontsize=7, ncol=2, frameon=False)
                status["figure_7"]["files"].append(_save(fig, folder, f"figure_7_{name}"))
        candidates_path = root / "candidates.csv"
        selections_path = root / "selections.csv"
        if candidates_path.exists() and selections_path.exists():
            candidates, selections = pd.read_csv(candidates_path), pd.read_csv(selections_path)
            xy = (
                ("projection_x", "projection_y")
                if "projection_x" in candidates
                else ("feature_0", "feature_1")
            )
            if set(xy).issubset(candidates.columns) and "gene_id" in candidates:
                for name, part in _strata(candidates):
                    query = selections
                    for key in ("scenario", "dataset"):
                        if key in part and key in query:
                            query = query[query[key].astype(str) == str(part.iloc[0][key])]
                    if query.empty:
                        continue
                    # One explicitly identified trajectory, not a union across seeds.
                    seed = query.seed.min() if "seed" in query else None
                    if seed is not None:
                        query = query[query.seed == seed]
                    if "policy" in query:
                        policy = (
                            "calder"
                            if (query.policy == "calder").any()
                            else str(query.policy.iloc[0])
                        )
                        query = query[query.policy == policy]
                    else:
                        policy = "saved"
                    if "threshold_quantile" in query:
                        quantile = sorted(query.threshold_quantile.unique())[
                            len(query.threshold_quantile.unique()) // 2
                        ]
                        query = query[query.threshold_quantile == quantile]
                    fig, ax = plt.subplots(figsize=(6.7, 5.2))
                    ax.scatter(part[xy[0]], part[xy[1]], s=16, c="0.8", label="Candidate features")
                    measured = part[part.gene_id.isin(query.gene_id)]
                    ax.scatter(
                        measured[xy[0]],
                        measured[xy[1]],
                        s=25,
                        facecolors="none",
                        edgecolors="tab:blue",
                        label="Measured in replay",
                    )
                    hit_ids = (
                        query.loc[query.is_hit.astype(bool), "gene_id"] if "is_hit" in query else []
                    )
                    hits = part[part.gene_id.isin(hit_ids)]
                    ax.scatter(
                        hits[xy[0]],
                        hits[xy[1]],
                        s=38,
                        marker="*",
                        c="tab:red",
                        label="Confirmed measured hit",
                    )
                    ax.set(
                        xlabel=xy[0].replace("_", " "),
                        ylabel=xy[1].replace("_", " "),
                        title=f"{name} | {policy}, seed={seed}\nProjected geometry does not establish causal biological relationships",
                    )
                    ax.legend(fontsize=8, frameon=False)
                    status["figure_8"]["files"].append(_save(fig, folder, f"figure_8_{name}"))
    for item in status.values():
        if item["files"]:
            item["status"] = "generated"
            item.pop("reason", None)
    (root / "figure_status.json").write_text(json.dumps(status, indent=2) + "\n")
    return status
