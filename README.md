# CALDER

**Calibrated Active Learning for Diverse and Efficient Biological Research**

CALDER is a reproducible Python research system for a finite-pool experimental-design question: **can probability calibration and diversity-aware, budget-constrained selection discover more high-effect perturbations, or cover more predefined groups, than standard selection rules?** Hit discovery and group coverage are separate primary outcomes. Calibration and robustness are evaluation properties; there is no combined headline performance score.

The repository implements synthetic experiments with known ground truth, a Bayesian surrogate, eight acquisition policies, leakage tests, backed AnnData preparation, paired evaluation, ablations, and figures generated from saved outputs. CALDER is a proposed combination to investigate, with no assumed superiority or new theorem.

No laboratory experiments are performed. Biological experiments, when supplied, are retrospective replays of published measurements. Synthetic results are never presented as biological validation.

## Installation

Use Python **3.11–3.13** and [uv](https://docs.astral.sh/uv/). No GPU, paid API, or cloud account is required. Run commands from the repository root:

```bash
uv sync --extra dev --frozen
uv run --frozen --extra dev calder doctor
uv run --frozen --extra dev pytest
uv run --frozen --extra dev ruff check .
uv run --frozen --extra dev ruff format --check .
```

`uv.lock` pins the resolved dependency graph; `pyproject.toml` bounds supported dependencies. Keep the lock file when sharing a run. For an activated environment, `source .venv/bin/activate` makes `calder` available directly. Optional saved-result exploration uses `uv sync --extra dev --extra demo --frozen`.

## Run a synthetic experiment

A small CPU configuration exercises the complete pipeline:

```bash
uv run --frozen --extra dev calder synthetic \
  --config configs/synthetic_smoke.yaml --run-dir results/smoke
uv run --frozen --extra dev calder report --run-dir results/smoke
```

The smoke configuration uses 60 candidates, two scenarios, two seeds and one threshold. Its calibration reserve is intentionally too small for the default minimum: raw-probability fallback warnings are expected. It is a software check, not an adequately replicated scientific comparison.

The full preregistered synthetic configuration runs four scenarios, 10 seeds, three thresholds and eight policies, plus an evaluation-only full-information reference:

```bash
uv run --frozen --extra dev calder synthetic --config configs/synthetic.yaml
```

The default pool has **250 candidates, 24 features, five groups and 128 random Fourier features**. Scenarios are A: smooth response, B: separated high-effect modes, C: heteroscedastic/noisy observations, and D: corrupted groups with distribution shift. Features are correlated within groups; latent effects, noisy measurements and predictions remain distinct. Outlier rate, overlap, corruption, noise and simulated heterogeneous costs are configurable.

The synthetic pilot uses 20% of the pool so that its 35% calibration reserve can support calibration. The real-data protocol uses a 5% pilot. Subsequent cost budgets are approximately 5%, 10% and 20% of pool size beyond the pilot, with a maximum batch size of 10. Fractional nominal budgets can leave an unspendable remainder when costs are one. A run directory is never overwritten; omit `--run-dir` for a timestamped directory, or choose a new empty destination for each replay.

## Data availability and biological workflow

The documented primary and secondary datasets are scPertEval's **Replogle 2022 essential-gene panels**, `replogle22k562` and `replogle22rpe1`. They are distinct from the larger genome-wide K562 screen. See the [official dataset documentation](https://scperteval.readthedocs.io/en/latest/user-guide/datasets.html).

On 3 November 2026, this execution environment could read the documentation, but direct object and public metadata requests for both H5AD sources returned HTTP 400 after an earlier sandbox DNS failure. This is an access limitation of the environment, not evidence that the public datasets do not exist. The record is in `data/manifests/source_verification.json`. **The K562 and RPE1 scientific benchmarks remain pending verified biological input.** Their outcomes and context comparison must not be inferred from synthetic results.

Inspect and verify before downloading:

```bash
uv run --frozen --extra dev calder data inspect --dataset replogle22k562
uv run --frozen --extra dev calder data download --dataset replogle22k562 --max-gb 2
uv run --frozen --extra dev calder data inspect --dataset replogle22rpe1
uv run --frozen --extra dev calder data download --dataset replogle22rpe1 --max-gb 2
```

Downloads stream to disk, enforce a size limit and record source/checksum information. Increase `--max-gb` only after checking the verified object size and available space. The source pattern is `https://storage.googleapis.com/scperteval/processed/<dataset>_processed_complete.h5ad`. An unavailable source produces an actionable error; there is no automatic synthetic substitute.

A verified local file is equally supported:

```bash
uv run --frozen --extra dev calder data inspect \
  --dataset replogle22k562 --path /absolute/path/replogle22k562.h5ad
uv run --frozen --extra dev calder data prepare \
  --dataset replogle22k562 --path /absolute/path/replogle22k562.h5ad \
  --config configs/replogle_k562.yaml
uv run --frozen --extra dev calder benchmark --config configs/replogle_k562.yaml

uv run --frozen --extra dev calder data prepare \
  --dataset replogle22rpe1 --path /absolute/path/replogle22rpe1.h5ad \
  --config configs/replogle_rpe1.yaml
uv run --frozen --extra dev calder benchmark --config configs/replogle_rpe1.yaml
```

Before preparation, edit a copy of the relevant config if your schema differs. `dataset.perturbation_key`, `control_label`, `normalization`, `local_path` and `prepared_path` declare the input convention and destinations. The defaults expect `obs['perturbation']`, the control label `control`, unique gene identifiers in `var_names`, and log-normalized `X`. Inspect these assumptions on the actual file. Gene symbols and Ensembl IDs are not silently treated as interchangeable.

`normalization: log_normalized` preserves `X`. For validated raw counts, set `normalization: counts`; the pipeline applies library-size normalization to 10,000 and `log1p` once. Raw-count cell QC is skipped and recorded as unavailable for already normalized input. Minimum perturbation size is 12 cells; the count-scale minimum is 200 expressed genes per cell. Response genes pass a detection filter in controls and are ranked by **control-only expression variance**, up to 1,024 genes. This is a documented variance-selection procedure, not a claim to reproduce a particular scanpy HVG method.

The loader uses backed AnnData and chunked/sparse access. It excludes ambiguous, combinatorial, unmapped and under-supported perturbations, records reasons, and materializes only bounded subsets required for statistics, bootstrap means and control embeddings. Keep raw/processed data outside Git; large files are ignored.

### Effect definition and uncertainty

For perturbation gene `g`, the response is

```text
delta_g[j] = mean_expression_g[j] - mean_expression_control[j]
effect_g = sqrt(sum_j weight[j] * delta_g[j]^2)
a_g = log(effect_g + epsilon)
```

Weights are stabilized, clipped inverse control variances normalized to sum to one. Set `dataset.effect_weighting: uniform` in a separate preparation config for an effect-definition sensitivity check. The score measures within-dataset transcriptional response. It is not a clinical phenotype, therapeutic benefit, or proof of a direct regulatory mechanism.

Seeded cell-level multinomial bootstrap resamples both perturbation and control cells. The prepared table stores effect magnitude, effect standard error, log-effect sampling variance, cell count, quality status and provenance. Fewer than 30 cells triggers a warning. A common control reference induces dependence between perturbation scores; the diagonal-noise surrogate does not model that dependence.

Each prepared dataset contains `pool.npz`, separately held `outcomes.npz`, `perturbations.csv` and `manifest.json`. The manifest records source URL, local checksum, exclusions, normalization, control-only transformations, response genes and weights, seeds, software processing version, feature/group semantics and artifact hashes. The real-data latent response is unknown: public measured effects are the offline reference, not error-free biological truth.

### Features and fixed groups

Control features comprise baseline mean, log variance, detection frequency and seeded truncated-SVD embeddings of control expression with genes as observations. Scaling and clustering use those permitted representations. The fallback groups are **artificial feature clusters**, not biological pathways.

The optional network API accepts a local, versioned human STRING TSV. Use explicit `gene1/gene2/score` columns, or STRING `protein1/protein2/combined_score` with a separate unambiguous `protein_id/gene_id` mapping. Supply the release version, preserve the input files and report unmapped candidates. Network features include degree, weighted degree and a bounded adjacency embedding; fixed-seed Louvain yields **network modules**. These are prior functional-association groups, not ground-truth causal pathways. The preparation API supports control, network, combined and uninformative feature modes. No cross-species bulk download is required.

For CLI preparation, a copied dataset config can declare:

```yaml
dataset:
  name: replogle22k562
  feature_mode: combined
  network_path: /absolute/path/human_candidate_edges.tsv
  network_version: "12.5"
  network_mapping: /absolute/path/protein_to_gene.tsv
```

Keep the other dataset settings from the supplied config. Omit `network_mapping` only when edge columns already contain exactly matching gene identifiers. Use a separate `prepared_path` and preparation `--output-dir` for each feature/effect variant.

## Frozen K562–RPE1 validation

After both datasets have been prepared, use the same method settings in the two configs:

```bash
uv run --frozen --extra dev calder validate-contexts \
  --k562-config configs/replogle_k562.yaml \
  --rpe1-config configs/replogle_rpe1.yaml \
  --output-dir results/context_validation
```

This command checks identical surrogate, acquisition and protocol settings and records a frozen-method hash. K562 is run first, then the same sequential procedure is rerun in RPE1 using a new measured pilot and subsequently revealed RPE1 observations. It builds the combined report in `results/context_validation/comparison`. This is **context re-evaluation**, not transfer of fitted K562 coefficients. Direct predictive transfer is not implemented. Do not select settings using final RPE1 outcomes. The context-comparison report requires genuine saved outputs from both contexts; missing outputs are reported as unavailable.

## Surrogate, calibration and selection

The principal surrogate is Gaussian-prior Bayesian regression over seeded random Fourier features. An intercept is included with the same prior precision. For design matrix `Phi`, observation sampling variances `s_i²` and configured variance floor `noise_floor`,

```text
v_i = max(s_i² + noise_floor, 1e-10)
A = prior_precision * I + Phi.T @ diag(1 / v_i) @ Phi
posterior_mean = solve(A, Phi.T @ (y / v))
latent_mean(x) = phi(x).T @ posterior_mean
epistemic_variance(x) = phi(x).T @ solve(A, phi(x))
```

Cholesky/triangular solves avoid explicit inversion. The runner converts the variance floor to the model's `observation_noise` standard-deviation parameter, so the variance is added once. Predictive variance returned by the surrogate is **epistemic variance of the latent log effect**; observation noise is not added to latent hit probability. The raw probability is `NormalCDF((mean - threshold) / sqrt(epistemic_variance))`, with explicit degenerate-variance handling.

Available models are global mean, deterministic ridge, fixed-prior Bayesian ridge, Bayesian random Fourier features, exact small-data GP, and a five-member bootstrap MLP ensemble. Ensemble disagreement/function draws are heuristic uncertainty, not an exact Bayesian posterior; the implementation uses CPU scikit-learn MLPs. Deterministic ridge returns zero epistemic variance. Main comparisons keep the same surrogate and parameters across all acquisition rules.

A uniformly random pilot is split into surrogate-training and calibration portions. The threshold is a predeclared 70th, 80th or 90th percentile **of pilot-training log responses only**; strict `>` defines a hit. Threshold and pilot IDs are identical across paired policies. Calibration reserve rows never enter surrogate training, and the surrogate is refitted using initial training rows plus newly acquired rows.

After each refit, isotonic calibration is refitted against the same fixed pilot reserve, with no actively selected labels added to the calibrator. Too few examples, a single label class or constant predictions triggers an explicit fallback to raw probabilities. Isotonic calibration targets **measured-hit labels**, whereas raw epistemic probabilities target latent effects. The distinction matters under measurement noise: the mapping does not establish latent-effect calibration. Adaptive selection also prevents an exchangeability guarantee.

The policies are:

| Policy | Decision rule |
| --- | --- |
| `random` | Uniform random candidate ordering with feasibility checks |
| `greedy_mean` | Largest predicted mean per unit cost |
| `uncertainty` | Largest epistemic standard deviation per unit cost |
| `ucb` | Mean plus configured multiple of epistemic standard deviation |
| `thompson` | One joint latent posterior function draw, ranked per cost |
| `probability_of_hit` | Raw probability of threshold exceedance per cost |
| `poh_diversity` | Raw hit probability plus expected new-group coverage |
| `calder` | Calibrated hit probability, expected new-group coverage, exploration and pair similarity penalty |

Thompson uses exact joint weight draws for Bayesian ridge/RFF, joint latent GP draws, or one bootstrap member for the ensemble. The general policy interface permits a documented independent-marginal fallback if no joint draw is supplied. PoH is an **independently implemented conceptual baseline**, not a claimed reproduction of the 2026 paper's complete method or benchmark protocol.

For batch `B`, CALDER greedily adds the greatest marginal utility per cost:

```text
sum_g p_hit(g)
+ diversity_weight * sum_{undiscovered k} [1 - product_g (1 - membership(g,k)*p_hit(g))]
+ uncertainty_weight * sum_g normalized_epistemic_std(g)
- similarity_weight * sum_{unordered pairs g,h} max(0, cosine(z_g,z_h))
```

Epistemic standard deviation is normalized by its maximum over the currently available pool. The group term is a working independent-hit approximation; correlated posterior responses and overlapping groups violate a literal independence interpretation. Previously hit groups receive no new-group reward. The pair penalty can make utility nonmonotone, so greedy selection has no asserted optimality guarantee and can stop before filling a batch if all feasible additions reduce utility. Ties follow fixed candidate order. Every batch enforces uniqueness, size and remaining cost. Biological costs are one; heterogeneous-cost simulations are explicitly labelled simulated.

## Leakage prevention and evaluation

`CandidatePool` contains features, costs and fixed memberships, with no responses. `OfflineOracle` provides selected queries, never a response-table property. Policy methods receive only unmeasured candidate metadata, predictive distributions and observed history. Full-label evaluation occurs after an entire trajectory has finished. Feature selection, scaling, SVD and group construction do not use unqueried intervention responses.

Tests deliberately alter hidden outcomes, compare subsequent selections, check reserve/training disjointness and reject duplicate queries. Python interfaces are research safeguards, not a security boundary against malicious introspection.

Primary metrics are cumulative confirmed measured hits, hit precision, recall relative to the offline measured pool, distinct predefined groups with confirmed hits, and discovery-versus-cost AUC. Counts include the common pilot; active-only counts are saved separately. Synthetic latent hits and latent group coverage are separate ground-truth diagnostics. AUC is the left-step integral `sum(diff(additional_cost) * preceding_hit_count)` over saved batch checkpoints: hits count only after their batch is measured. Its units are hit × cost; normalization by the realized additional-cost span gives the mean number of available hits during acquisition. Pilot acquisition cost is excluded from that integration, while pilot hits form the starting count. Consult run outputs for realized costs and unspent budget.

The offline full-information ranking is evaluation-only. With equal costs it is a hit-count upper-bound reference for the measured-label objective; it is not deployable and **does not upper-bound group coverage**. The heterogeneous-cost ranking is a reference, not a knapsack-optimal bound.

Secondary diagnostics include prequential Brier score, fixed-bin ECE, mean absolute prediction error, uncertainty–error association, interval coverage, runtime, memory where recorded and seed variability. All prequential predictions are saved **before** corresponding outcomes are revealed. Diagnostics on actively acquired examples are selection biased and do not establish calibration on the full pool. Raw and calibrated predictions can have different biological/statistical targets; results must state the target.

Bootstrap intervals resample whole paired seed runs within each dataset/scenario/threshold. Seed-level rows remain available. Repeated seeds quantify algorithmic variation on a fixed pool, **not independent biological replications**. Multiple descriptive contrasts do not justify an unadjusted universal superiority claim. Negative or mixed outcomes are retained.

## Ablations and sensitivity

```bash
uv run --frozen --extra dev calder ablate --config configs/ablations.yaml
uv run --frozen --extra dev calder ablate --config configs/ablations.yaml \
  --sensitivity --run-dir results/development_sensitivity
```

Component ablations remove calibration, diversity, uncertainty or similarity; randomize group assignments; or replace informative features. Randomized-group decisions are evaluated against the original fixed groups. The small declared development grid varies pilot size (10%, 20%, 30%), batch size (5, 10, 20), simulated noise (0.05, 0.2, 0.5) and surrogate (Bayesian ridge, RFF, ensemble). Threshold sensitivity is already crossed through the declared quantiles. All declared variants are retained; there is no automatic winner-driven tuning. Use synthetic/development runs to choose and freeze a method before final RPE1 evaluation.

Network/control/combined feature comparisons require separately prepared, provenance-matched biological pools. Alternative effect definitions likewise require separate preparation. These are not fabricated by the synthetic component-ablation command. Optional E-distance is a bounded, independently implemented empirical score available through `calder.data.pseudobulk.energy_distance`; it does not replace the primary effect definition.

## Saved outputs and figures

Every run saves its resolved `config.json`, configuration hash, environment/dependency versions, input hashes and Git commit/dirty-state information when Git is available. To replay, reuse the exact saved config, lock file, input data and code revision. Timestamps and runtimes are expected to differ; scientific trajectories should match for a fixed environment and seed.

| Artifact | Meaning |
| --- | --- |
| `experiment_manifest.json` | Run status, provenance, config/input hashes and protocol descriptions |
| `rounds.csv` | Discovery, coverage and realized budget at each logged round |
| `selections.csv` | Pilot/active selections, outcomes, costs and decision components |
| `prequential.csv` | Selected-candidate predictions recorded before reveal; evaluation diagnostics appended afterwards |
| `seed_level_results.csv` | Every seed/scenario/threshold/policy endpoint |
| `results.csv` | Aggregated endpoint estimates and bootstrap intervals |
| `paired_differences.csv` | Paired policy contrasts and intervals |
| `calibration.csv` | Raw/calibrated diagnostics with explicit measured/latent targets |
| `ablations.csv` | Component-ablation endpoints, or a header-only table in ordinary runs |
| `candidates.csv` | Outcome-free coordinates for saved-result projection figures |
| `trajectories/*/events.jsonl` | Query order, training/reserve IDs and audited reveal events |
| `trajectories/*/predictions_*.npz` | Full available-pool prediction snapshots written before each reveal |
| `summary.md`, `figures/`, `figure_status.json` | Regenerable report, plots and explicit figure availability |

The eight requested figure families cover discovery, group coverage, reliability, discovery–coverage trade-offs, ablations, K562–RPE1 comparison, uncertainty versus error, and candidate geometry. A report generates only figures supported by saved inputs. An ordinary synthetic run cannot produce a biological context comparison; an ordinary non-ablation run cannot provide ablation evidence. Missing figures are recorded as unavailable, not filled with example values. Projected feature distances do not establish causal biological relationships.

The exploratory notebook and optional demo read saved outputs:

```bash
uv run --frozen --extra dev --extra demo calder demo --run-dir results/smoke
```

## Repository map and quality checks

`src/calder/data` handles provenance, validation and sparse preparation; `models` contains surrogate inference/calibration; `acquisition` contains outcome-free decisions; `experiments` contains synthetic ground truth, the restricted oracle and paired runner; `evaluation` handles metrics, bootstrap statistics and plots. Scripts mirror the CLI. `notebooks/exploratory_analysis.ipynb` is a read-only saved-result explorer. The scientific narrative and bibliography are in `reports/`.

Unit and integration tests cover analytical Gaussian posteriors, epistemic variance, joint draws, degenerate matrices, calibration fallbacks, costs, empty pools, group coverage, sparse/backed data, bootstrap aggregation, hidden-label isolation, serialization, audit logs and reproducibility. Run `pytest` before changing a milestone status. `calder doctor` checks Python/dependencies, output access, disk guidance, configured prepared data and a minimal synthetic replay; GPU support is optional.

## Scientific interpretation and sources

A large expression-response score does not establish therapeutic usefulness. Perturbation-induced changes do not automatically identify direct causal regulatory edges. STRING communities are prior association groups. Retrospective replays omit real prospective assay scheduling and many forms of biological/technical shift. Calibration can improve probability scores while leaving decisions unchanged; diversity can increase coverage while reducing hit discovery; PoH or another baseline can win. Those outcomes are scientifically informative.

The verified reading record and any source-access limits are in `reports/bibliography.md`. Core references are:

- Rubbi et al. (2026), [Many Needles in a Haystack: Active Hit Discovery for Perturbation Experiments](https://arxiv.org/abs/2605.10196).
- [Deep learning perturbation models can outperform baselines on calibrated metrics](https://www.nature.com/articles/s41587-026-03307-w), Nature Biotechnology (2026).
- [scPerturb: harmonized single-cell perturbation data](https://www.nature.com/articles/s41592-023-02144-y), Nature Methods.
- [scPertEval dataset documentation](https://scperteval.readthedocs.io/en/latest/user-guide/datasets.html).
- [STRING functional association networks](https://string-db.org/).

See `CITATION.cff` for software citation metadata and `LICENSE` for the MIT license.
