# CALDER: Uncertainty-Calibrated, Diversity-Aware Active Learning for Budgeted Genetic Perturbation Discovery

## Abstract

We implemented a leakage-isolated, CPU-based system for budgeted perturbation selection using Bayesian random-feature regression, pilot-reserved probability calibration and batch diversity. Controlled evaluation comprised 1,080 primary synthetic trajectories and 840 component-ablation trajectories. We also completed 270 retrospective trajectories in each of the Replogle K562 essential and RPE1 panels, using the checksum-verified scPerturb release and a frozen method. At the 80th-percentile pilot threshold, CALDER found fewer synthetic hits than Probability-of-Hit in all four landscapes. In biological replay, CALDER minus Probability-of-Hit was −14.3 confirmed hits [95% paired bootstrap CI −18.0, −10.6] in K562 and −20.1 [−28.0, −10.3] in RPE1; artificial control-feature cluster coverage changed by +0.3 [0.1, 0.6] and +0.3 [0.0, 0.6] respectively. These clusters are diversity proxies, not validated biological pathways. Calibration improved acquired-label Brier score in K562 despite poorer hit discovery, illustrating that probability reliability and decision efficiency can diverge. RPE1 calibration changed little and predictive intervals undercovered in both contexts. Results describe ten paired algorithmic seeds per fixed context, not biological replications or laboratory experiments performed by CALDER. The evidence supports retaining negative outcomes and reporting hit discovery, proxy-group coverage and calibration separately; it does not demonstrate an overall CALDER advantage.

## 1. Introduction

The experimental-selection problem is to choose a finite sequence of single-gene perturbations under a fixed experimental budget. The primary outcomes are the number of confirmed high-effect perturbations and the number of predefined groups containing confirmed hits. Prediction quality and calibration are evaluation properties rather than components of a composite headline score.

CALDER studies whether calibrated threshold-exceedance probabilities, group diversity and an optional exploration bonus change the efficiency of sequential selection. It is an empirical extension to existing active-discovery ideas, not an asserted novel theorem. No real laboratory experiments are performed by this software.

## 2. Related work

Probability-of-Hit motivates selecting interventions by posterior threshold-exceedance probability [1]. The baseline here is independently implemented, sharing CALDER's surrogate so that the primary comparison isolates the decision rule. It does not reproduce the full published protocol.

Perturbation benchmark metrics can be insensitive to useful predictive signal [2]. Metric calibration in that work is distinct from probability calibration here. Harmonized public perturbation resources enable auditable retrospective analysis [3,4]. STRING provides prior associations useful for fixed network modules [5]; associations are not causal regulatory edges. Replogle and colleagues performed the original Perturb-seq experiments [6]. See [bibliography](bibliography.md) for source metadata and verification limits.

## 3. Research questions

1. At equal measured experimental cost, does CALDER discover more confirmed hits than independently implemented Probability-of-Hit using the same surrogate?
2. Does expected group coverage increase final observed group coverage, and what is the associated change in hit count?
3. Does calibrating probabilities reduce prequential Brier score or ECE, and does any improvement alter selection quality?
4. How do these differences vary under heteroscedastic noise, representation misspecification, unreliable group labels and cell-context shift?

The null includes no improvement, and negative outcomes remain in all output tables. A monotone probability calibration can leave single-candidate rankings unchanged. In a diversified batch objective it can change the relative scale and saturation of coverage utility. Calibration quality and decision quality must therefore be evaluated separately.

## 4. Methodology

### Response definition

For perturbation g, subtract the untreated-control mean from its cell-mean expression vector. Define `effect(g) = sqrt(sum_j w_j delta_gj²)`. Weights are nonnegative stabilized, clipped inverse control variances, normalized to sum to one. The implementation records the exact epsilon and clipping configuration. This is a within-dataset transcriptional-response magnitude, not a validated clinical phenotype.

The model target is `a_g = log(effect(g) + epsilon)` with the threshold on the same log scale. Cell-level seeded bootstrap estimates sampling uncertainty. Both perturbation cells and controls should be resampled when estimating uncertainty relative to a sampled control population. Shared controls induce correlated response-estimation errors; the surrogate's diagonal observation-noise approximation does not model that correlation. Bootstrap standard errors do not eliminate biological, batch or guide-efficiency confounding.

### Permitted representations

Feature selection and transformations use untreated-control cells or fixed external annotations only. Control-based means, variances, detection rates and gene embeddings are available before interventions are queried. Fixed feature clusters are artificial clusters; versioned mapped STRING communities are network modules. Missing identifiers and exclusions must be reported. None of these labels is an experimentally validated pathway assignment.

### Surrogate

Bayesian regression over seeded random Fourier features uses Gaussian prior `w ~ N(0, alpha^-1 I)` and heteroscedastic diagonal observation variances. For feature matrix Phi and precision matrix W, posterior precision is `A = alpha I + Phi.T W Phi` and posterior mean solves `A m = Phi.T W y`. Cholesky solves replace explicit matrix inversion. Epistemic prediction variance is `phi_g.T A^-1 phi_g`; observation variance is accounted for separately.

The latent hit probability is `P(a_g >= tau | observed history) = 1 - Phi_Normal((tau-mu_g)/sigma_epistemic_g)`. Near-zero variance has a safe deterministic limit. Ridge, Bayesian linear, mean, optional GP and optional ensemble variants expose the same prediction interface. Predictive accuracy alone need not imply acquisition efficiency.

For measured-response predictive intervals, variance is the sum of epistemic variance, estimated cell-sampling variance and the configured observation-model discrepancy floor. These intervals are labelled `measured_predictive_model_assumed`. The floor is a likelihood assumption, and this interval construction is heuristic for ridge/ensemble uncertainty. Latent-effect intervals use epistemic variance alone. Older saved trajectories lacking a separate discrepancy column restore it from the resolved configuration in memory; their original prediction CSVs are preserved.

### Calibration

A uniformly selected pilot is split into disjoint surrogate-training and calibration reservations. The calibrator is fitted against held-out pilot hit labels, while subsequent surrogate fits exclude calibration-reserved candidates. Tiny or single-class calibration sets trigger a documented raw-probability fallback. Reusing a fixed pilot calibration set after surrogate updates is a practical heuristic rather than a guarantee of validity under a changing score distribution.

Calibration against noisy pilot hit indicators estimates a measured-label correction. It is not guaranteed to calibrate the probability of exceeding a latent biological threshold. Synthetic data separately evaluate latent and measured labels to expose this distinction.

### Batch objective and optimization

CALDER maximizes a greedy approximation to `sum p_hit + lambda expected_new_group_coverage + eta uncertainty_bonus - rho within_batch_similarity`. Previously hit-containing groups are excluded from new-coverage utility. For each remaining group k, the approximation is `1 - product_g (1 - membership_gk p_hit_g)`. This assumes conditional independent hit indicators and treats membership weights as probabilities in that term; correlated biological effects can violate the approximation.

Candidates are chosen by marginal utility per unit cost when costs differ. Batch size, remaining budget, no replacement and deterministic tie-breaking are enforced. Similarity penalties mean the complete objective need not be monotone; greedy selection has no asserted global-optimality guarantee. Biological primary runs use equal costs; heterogeneous costs are explicitly simulated.

### Leakage isolation

Candidate metadata contain no response labels. Acquisition receives only candidate features, predictive distributions, observed history, fixed group assignments and costs. The oracle reveals queried observations only. Full-pool labels are used by a separate evaluator after decisions. Tests perturb hidden outcomes and verify identical next decisions, reserve disjoint calibration data, reject repeated selection and replay seeded trajectories.

## 5. Experimental setup

Synthetic scenarios cover a smooth landscape, separated hit modes, heteroscedastic observations, and inaccurate groups with distribution shift. Latent effects, noisy measurements and predictions remain distinct. Repeat over predeclared seeds with paired pilots for all policies. The synthetic configuration is the laptop-sized end-to-end entry point.

Biological runs use the Replogle essential-gene K562 panel followed by independent RPE1 reruns with frozen method hyperparameters. Uniform pilots supply preregistered percentile thresholds; the hidden complete outcome pool never sets the threshold. Report approximately 5%, 10% and 20% additional-budget checkpoints where feasible. Distinguish initial pilot expense from additional acquisition cost.

Baseline policies are random, predicted mean, epistemic uncertainty, UCB, Thompson sampling, Probability-of-Hit, Probability-of-Hit with diversity, and CALDER. The full-information offline ranking is only an upper-bound reference and cannot acquire observations for a deployable policy.

The experiment manifest records resolved configuration, config hash, source and processed-data provenance, seeds, software versions, Git revision and working-tree state. Repeated seeds on one fixed dataset measure algorithmic variation, not independent biological replication.

### Completed synthetic runs and provenance

The primary run is [`results/reference_synthetic`](../results/reference_synthetic/summary.md): 250 candidates, 24 feature dimensions and five predefined synthetic groups per scenario; data seed 317; acquisition seeds 0–9. The shared pilot contains 50 measured candidates (33 surrogate-training, 17 calibration-reserved). Batch size is at most ten, and additional unit-cost checkpoints are 12.5, 25 and 50 available budget units; all final trajectories measured 50 additional candidates, so final metrics include 100 total experiments including the pilot. Thresholds are the 70th, 80th and 90th percentiles of pilot-training responses. The q=0.8 summaries below were chosen as the central predeclared threshold, not because of their results; the other thresholds remain in the saved outputs.

The main surrogate has 128 random Fourier features, prior precision 1.0, length scale 2.0 and observation-variance floor 0.04. CALDER weights are diversity 1.0, uncertainty 0.1 and similarity 0.1. These are reference settings, not a demonstrated optimum. Scenario A has a smooth linear log-response; B has separated nonlinear hit modes; C adds heteroscedastic observations; D combines context shift and corrupted groups. One seeded pool exists for each scenario. The four scenarios are not treated as four independent replications of one population effect.

The primary configuration hash is `a76b99c6129af0d9adc3ed4cb1beaf456ba1beb0fa75f28503d87191a1ec2a52`; the ablation configuration hash is `ffe64235d97614f80d76e8a075a8d1ea29c7bccd5d185ddf3747e27cba281a81`. Both manifests record Git `15de4bebfb212eca645a74806bafbee35beb3a5b` with `git_dirty=true`. The runner creates the untracked result directory before recording this flag; the flag therefore does not by itself establish a tracked-source modification. The implementation session reported unchanged tracked execution source at run startup. The original flag is preserved, and the delivered source, resolved configurations and manifests should be retained together. The primary `seed_level_results.csv` SHA-256 is `8012c6cb4afda2aa45a27c64dfd965179621b7665e85233fe2d486671d0651f8`, and the ablation endpoint-table SHA-256 is `6a43bec664085c5b587ee5ee018c26814bbc696e9d73db3bfda4bdca86aa2821`. These four key saved files (configuration, seed endpoints, paired differences and calibration diagnostics) were checked against each run manifest. There were no repeated candidate indices within any saved prequential trajectory.

Primary acquisition/evaluation runtime was 29.67 seconds; the ablation run was 27.14 seconds on the recorded CPU environment, excluding subsequent figure/report generation. Process lifetime resident high-water marks were 492.3 MB and 443.9 MB respectively; these include imports and are not incremental per-trajectory allocations.

### Completed biological preparation and frozen context validation

The exact source is scPerturb Zenodo release 1.4, [record 13350497](https://zenodo.org/records/13350497), not the scPertEval processed files. The K562 **essential** file and RPE1 file were downloaded with verified publisher MD5 checksums. Backed/chunked count-scale processing applies library-size normalization to 10,000 and log1p once. Newly fitted expression filtering, response-gene selection, statistical features, SVD and clustering use untreated controls only.

| Context | Source cells / genes | Retained controls | Retained intervention cells | Excluded cells | Candidate genes | Low-cell candidates (12–29 cells) |
| --- | --- | --- | --- | --- | --- | --- |
| K562 essential | 310,385 / 8,563 | 10,691 | 274,151 | 25,543 | 1,855 | 66 |
| RPE1 | 247,914 / 8,749 | 11,485 | 201,802 | 34,627 | 2,045 | 290 |

Each context uses 1,024 response genes chosen by control-only variance after control detection filtering, 16 gene-embedding components from a bounded 2,048-control-cell subset, and five context-specific artificial control-feature clusters. No STRING network was supplied, so coverage findings concern these proxies rather than independently validated biological modules. Single-gene candidates require at least 12 retained cells and an identifier represented in the released gene universe. Fifty seeded bootstrap replicates resample perturbation and control cells; 12–29-cell candidates remain with explicit uncertainty warnings.

The biological pilot is 5%: K562 has 93 pilot candidates (61 training, 32 calibration) and 371 additional experiments at the final 20% checkpoint; RPE1 has 103 pilot candidates (67 training, 36 calibration) and 409 additional experiments. Final precision therefore divides by 464 and 512 total experiments respectively, including pilots. Thresholds remain the predeclared 70th/80th/90th pilot-training percentiles. Each context has ten seeds, eight acquisition rules and an evaluation-only reference, for 270 trajectories.

The shared method was frozen before biological selection in [frozen_protocol.json](frozen_protocol.json), with verified method hash `d0556a6e2cc15b77a771a7fcc3ee111bf996aaf9477b6e5dc54511992ac3bfe7`. Model, acquisition and complete protocol sections match that frozen object in both context configurations. Both runs record Git `7e95d257bd7aaec8247dae6ff8f8d524fdb2c024`, `git_tracked_dirty=false`, and identical source-tree hash `a4055224cc2b528efbbedf5e60840df7ee89f4cf75338d07bbed4c9a83100abf`. Their `git_dirty=true` flag includes untracked generated outputs. Original experimental provenance is preserved separately from later report-generation provenance.

- K562: config `f22d66d39fdb00d01c7278aefdc90bc7bd0c5d229826df854e0b7227b475e1b1`; endpoint CSV SHA-256 `12d0a6ed01eeb972ecfd9f304640e2ca135ec0f81dd4bb41502af21efe496254`.
- RPE1: config `00f47d76e3508bbcd8248379fcb1b19b0128e1b981597bcaf959e7708ecc44ee`; endpoint CSV SHA-256 `053f4f785ad3ce1020971980686a8c96223526e2f912009f4323fe757d31ebc6`.

Detailed QC and processing provenance are in [K562 preparation](../data/manifests/replogle22k562_preparation.json) and [RPE1 preparation](../data/manifests/replogle22rpe1_preparation.json). Both runs' configuration, endpoint, paired-difference and calibration table hashes were checked against their manifests. This is independent sequential re-evaluation in RPE1 using newly revealed observations, **not** direct transfer of a fitted K562 model. No RPE1 outcome-label tuning was performed.

## 6. Results

### Hit discovery and predefined-group coverage

The table reports conditional means across the same ten paired seeds at the 80th-percentile pilot-training threshold, after 50 additional experiments. Confirmed hits use strictly greater-than-threshold measured log effects; groups are the five fixed **synthetic groups**, not biological pathways. Differences are CALDER minus Probability-of-Hit using the same random-feature surrogate, pilot split and threshold.

| Scenario | Confirmed hits: CALDER / PoH | CALDER − PoH hits [95% CI] | Groups: CALDER / PoH | CALDER − PoH groups [95% CI] |
| --- | --- | --- | --- | --- |
| A | 32.6 / 40.6 | -8.00 [-10.30, -5.80] | 3.1 / 3.1 | +0.00 [+0.00, +0.00] |
| B | 20.5 / 24.3 | -3.80 [-6.20, -1.40] | 5.0 / 4.7 | +0.30 [+0.10, +0.60] |
| C | 23.6 / 26.7 | -3.10 [-4.50, -1.80] | 4.9 / 5.0 | -0.10 [-0.30, +0.00] |
| D | 22.5 / 23.8 | -1.30 [-3.00, +0.50] | 5.0 / 5.0 | +0.00 [+0.00, +0.00] |

The primary hit-discovery hypothesis was not supported by these reference settings: CALDER had fewer measured hits in every scenario at q=0.8. In B, the additional 0.3 groups came with 3.8 fewer hits, quantifying a coverage–discovery trade-off rather than an overall win. In A, group coverage was identical despite eight fewer hits. C's groups were already nearly saturated, limiting any scope for a coverage gain. D's hit-count interval included zero. Bootstrap intervals resample whole paired seeds (2,000 resamples); they are descriptive conditional Monte Carlo intervals, with no biological-replication or multiplicity-adjusted superiority claim.

This conclusion is not limited to the displayed threshold: CALDER's mean measured-hit difference was negative in all twelve scenario/threshold strata. The full table, including every seed and the q=0.7/q=0.9 thresholds, remains in [paired differences](../results/reference_synthetic/paired_differences.csv) and [seed-level results](../results/reference_synthetic/seed_level_results.csv). No scenario-pooled effect estimate is reported.

Synthetic latent-hit means at q=0.8 were 31.9, 19.9, 21.4 and 18.3 for CALDER versus 40.4, 22.6, 24.0 and 20.4 for Probability-of-Hit in A–D respectively. Thus the measured-hit disadvantage was not reversed by consulting known latent truth after acquisition. The full-information reference ranks measured responses and bounds measured-hit count under these equal costs; it is not a latent-hit or group-coverage upper bound. For example, in noisy scenario C it found 57.6 measured hits but only 37.8 latent hits on average.

### Prequential probability diagnostics

Each cell below is the mean of ten seed-level scores on CALDER's own 50 adaptively acquired observations per seed at q=0.8. Raw and mapped probabilities are compared on the **same acquired examples** within each run; this is not a comparison across different policy-selected populations. Lower Brier score and ECE are better. Fixed-width ten-bin ECE is descriptive and depends on binning and finite sample size.

| Scenario | Target | Raw Brier | Mapped Brier | Raw ECE | Mapped ECE |
| --- | --- | --- | --- | --- | --- |
| A | measured | 0.175 | 0.292 | 0.171 | 0.332 |
| A | latent | 0.173 | 0.296 | 0.177 | 0.337 |
| B | measured | 0.219 | 0.309 | 0.259 | 0.313 |
| B | latent | 0.224 | 0.325 | 0.257 | 0.331 |
| C | measured | 0.198 | 0.289 | 0.173 | 0.261 |
| C | latent | 0.189 | 0.287 | 0.190 | 0.283 |
| D | measured | 0.222 | 0.326 | 0.225 | 0.338 |
| D | latent | 0.217 | 0.319 | 0.252 | 0.350 |

The monotone mapping worsened mean Brier score and ECE in every displayed scenario/target combination. It was fitted in every acquisition round for A, C and D, and in 90% of B's rounds; the remaining B cases used the documented raw-probability fallback. These results provide no evidence that this small-pilot calibration protocol improved probability reliability on the active selection distribution. The calibration reserve is only 17 observations, and changing surrogate fits can shift its probability-score distribution. These are plausible failure mechanisms to test, not established explanations.

Measured-label and latent-label diagnostics are kept separate because the map is fitted to noisy measured-hit labels whereas the raw Gaussian exceedance probability uses epistemic uncertainty for a latent log effect. Across q=0.8, latent epistemic 95% interval coverage was 0.966, 0.954, 0.990 and 0.958 in A–D. These descriptive coverages do not imply calibrated hit probabilities or exchangeability-based guarantees. Full metrics, bin counts and group-specific diagnostics are saved in [calibration.csv](../results/reference_synthetic/calibration.csv), [reliability.csv](../results/reference_synthetic/reliability.csv) and [grouped_calibration.csv](../results/reference_synthetic/grouped_calibration.csv).

### Probability-plateau diagnostic

A separate audit inspected all 280 saved CALDER q=0.8 pre-query snapshots (70 per scenario) without changing the policy or tuning its settings. Raw probabilities were distinct for every candidate in every snapshot. The mapped probabilities frequently flattened large parts of the pool:

| Scenario | Mean distinct raw values / pool size | Mean distinct mapped values / pool size | Mean largest equal-probability plateau / pool size | Entirely constant mapped pools |
| --- | --- | --- | --- | --- |
| A | 100% | 16.0% | 59.8% | 0/70 |
| B | 100% | 21.6% | 54.1% | 0/70 |
| C | 100% | 10.4% | 64.4% | 3/70 |
| D | 100% | 10.8% | 67.1% | 3/70 |

These are unweighted means across round snapshots; pools shrink as candidates are acquired. No snapshot was entirely zero or entirely one after mapping. Reconstructing every within-batch marginal objective from saved probabilities, permitted features and previously observed hit groups found a **unique** maximizing candidate in all 2,000 acquired decisions, with all saved selections matching that maximum. Thus deterministic row-order tie-breaking is not supported as an explanation for these trajectories under the nonzero uncertainty bonus. Probability flattening is observed; a resulting increase in the relative influence of uncertainty, diversity or similarity is a plausible mechanism, not a causal conclusion. Calibration removal tests change whole trajectories and do not isolate this mechanism by themselves.

The derived [audit table](synthetic_probability_plateau_audit.csv) and [input-provenance record](synthetic_probability_plateau_audit.json) are reproducible with `python scripts/audit_probability_plateaus.py --run-dir results/reference_synthetic --output reports/synthetic_probability_plateau_audit.csv`. The script checks regenerated feature hashes against the source run, audits saved decisions, and leaves the benchmark outputs unchanged.

### Retrospective K562 and RPE1 validation

The following conditional means use the central predeclared q=0.8 threshold at the final 20% additional-budget checkpoint. Precision and recall include the pilot; recall is relative to the context's complete offline **measured-hit** pool at each paired pilot-derived threshold.

| Context | Policy | Confirmed hits | Precision | Recall | Artificial clusters with hits | Total experiments |
| --- | --- | --- | --- | --- | --- | --- |
| K562 | CALDER | 95.9 | 20.67% | 25.52% | 4.5 | 464 |
| K562 | PoH | 110.2 | 23.75% | 29.42% | 4.2 | 464 |
| RPE1 | CALDER | 126.5 | 24.71% | 27.65% | 5.0 | 512 |
| RPE1 | PoH | 146.6 | 28.63% | 32.55% | 4.7 | 512 |

Paired differences are CALDER minus Probability-of-Hit; precision/recall differences are percentage points. The ten seeds are paired within each context, not across cell types.

| Context | Δ hits [95% CI] | Δ clusters [95% CI] | Δ precision, pp [95% CI] | Δ recall, pp [95% CI] |
| --- | --- | --- | --- | --- |
| K562 | -14.30 [-18.00, -10.60] | +0.30 [+0.10, +0.60] | -3.08 [-3.88, -2.28] | -3.89 [-4.95, -2.84] |
| RPE1 | -20.10 [-28.00, -10.30] | +0.30 [+0.00, +0.60] | -3.93 [-5.47, -2.01] | -4.90 [-7.15, -2.51] |

CALDER discovered fewer hits in both contexts, with hit-count intervals excluding zero, while K562 gained a small amount of artificial-cluster coverage. RPE1's corresponding coverage interval included zero. There is no demonstrated improvement in validated biological-group coverage because no independent network/pathway labels were used in these completed runs.

Absolute cross-context hit counts should not be read as a cell-type hit-rate comparison: RPE1 measured 512 total candidates and K562 464, their response scores are standardized within different datasets, their pilot-derived thresholds differ, and their feature clusters have different identities. Precision is a clearer budget-normalized descriptive comparison, but still concerns context-specific hit definitions. The [context comparison](../results/reference_contexts/summary.md) and [Figure 6](../results/reference_contexts/figures/figure_6_q-0.8_random_features.png) display actual saved counts and precision under a frozen procedure, with candidate and experiment counts labelled; they do not establish that one cell type is biologically easier to perturb.

### Biological calibration and model failure

These are means of seed-level diagnostics on CALDER's adaptively acquired measured outcomes, with both probability versions scored on the same acquired examples.

| Context | Raw Brier | Mapped Brier | Raw ECE | Mapped ECE | Model-assumed 95% interval coverage |
| --- | --- | --- | --- | --- | --- |
| K562 | 0.268 | 0.214 | 0.265 | 0.162 | 0.707 |
| RPE1 | 0.296 | 0.295 | 0.282 | 0.258 | 0.600 |

K562 shows a clear descriptive improvement in probability metrics despite fewer discovered hits, so improving prediction reliability did not guarantee an improved selection rule. RPE1's mean Brier change is very small (0.29584 to 0.29457), with a modest ECE reduction. Both contexts have substantial predictive-interval undercoverage even after adding sampling variance and the configured model-discrepancy floor: 0.707 and 0.600 against nominal 0.95. This suggests inadequacy of the assumed predictive distribution in these selected examples; it is not evidence of a guaranteed calibrated latent-effect model. True latent biological effects are unavailable in these real-data replays.

Source-dependent noise, guide efficacy, shared controls and context shift remain plausible contributors. RPE1 contains more low-cell candidates (290 versus 66) and a lower median retained cell count per candidate (70 versus 121), but these observations do not isolate a causal explanation. The [K562 uncertainty–error plot](../results/reference_k562/figures/figure_7_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png) and Figure 6 were visually checked for correct measured-target labels and legibility; dense scatter overlap should not be interpreted as a calibrated uncertainty–error relationship.

### Figure and source provenance

All figures derive from saved outputs. Figure 1A (q=0.8), Figure 3C (q=0.8) and Figure 5D (q=0.8) were visually inspected: labels, legends, axes and negative effects were legible. [Discovery curve](../results/reference_synthetic/figures/figure_1_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [reliability diagram](../results/reference_synthetic/figures/figure_3_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [ablation plot](../results/reference_ablations/figures/figure_5_dataset-synthetic__scenario-D__threshold_quantile-0.8.png). The reliability diagram pools acquired events descriptively; its individual bins are not independent validation samples. Figure 6 is now available from the completed K562/RPE1 saved outputs, with unequal experiment counts explicitly documented above.

The scPertEval GCS endpoints returned HTTP 400 in this execution environment. The verified alternative scPerturb release supplied `ReplogleWeissman2022_K562_essential.h5ad` (1,546,729,675 bytes) and `ReplogleWeissman2022_rpe1.h5ad` (1,236,886,900 bytes). Both downloads were checksum-verified and independently schema-validated as counts before normalization. They are not silently relabelled scPertEval processed files. Source URLs, checksums and access checks remain in [source_verification.json](../data/manifests/source_verification.json). Synthetic H5AD unit-test fixtures are distinct from the completed public-data biological replays.

## 7. Ablation studies

The separate [`reference_ablations`](../results/reference_ablations/summary.md) run contains 840 trajectories: full CALDER and six controlled variants over the same four pools, three thresholds and ten seeds. Full-CALDER means match the primary run. The following differences are **variant minus full CALDER** in confirmed hits at q=0.8; positive values therefore favour removing or replacing the named component.

| Variant minus full CALDER: hits [95% CI] | A | B | C | D |
| --- | --- | --- | --- | --- |
| Remove calibration | -0.50 [-2.00, +1.20] | +1.10 [-1.50, +3.70] | +2.10 [+0.60, +3.60] | +0.70 [-1.20, +2.10] |
| Remove diversity | +0.60 [-0.50, +1.90] | -0.30 [-0.90, +0.20] | +0.40 [+0.00, +1.20] | -0.80 [-2.20, +0.40] |
| Remove uncertainty bonus | -0.10 [-3.00, +2.80] | +0.20 [-0.80, +1.10] | +1.00 [-0.90, +2.80] | +0.40 [-0.60, +1.40] |
| Remove similarity penalty | +3.70 [+0.90, +6.30] | +1.30 [-0.40, +3.00] | +0.90 [-1.10, +2.70] | +0.30 [-1.60, +2.20] |
| Random groups | +0.20 [-1.20, +1.80] | +0.20 [-0.50, +1.10] | +0.10 [-0.90, +1.20] | -0.70 [-2.70, +1.00] |
| Uninformative features | -14.90 [-17.60, -12.40] | -1.40 [-2.80, +0.00] | +0.40 [-1.50, +2.40] | -1.10 [-3.80, +1.50] |

Removing the similarity penalty improved A by 3.7 hits [0.90, 6.30], with a group difference of +0.1 [0.00, 0.30]. The penalty therefore harmed discovery under this reference setting without a demonstrated coverage benefit there. Removing calibration improved noisy scenario C by 2.1 hits [0.60, 3.60]; its group difference was +0.1 [0.00, 0.30]. Calibration did not produce a consistent acquisition benefit in the other scenarios.

In separated-mode scenario B, removing diversity changed hit count by −0.3 [−0.90, 0.20] and group coverage by −0.4 [−0.80, 0.00]. Removing uncertainty changed hits by +0.2 [−0.80, 1.10] and groups by −0.3 [−0.80, 0.00]. These borderline/descriptive intervals do not establish a reliable component benefit, but they identify coverage as the relevant outcome rather than hit count. Under the fixed settings, randomized group labels also changed B's group coverage by −0.3 [−0.70, 0.00]. Uninformative representations strongly harmed the smooth A landscape (−14.9 hits), while effects in nonlinear/noisy scenarios were smaller and uncertain.

These are component tests at one fixed setting, not proof that an exhaustive search could not find a better setting. No hyperparameter was selected retrospectively to make CALDER win. The main component run does not include network-feature comparisons. Separate [exploratory model comparisons](../results/reference_models/summary.md) evaluated Bayesian ridge, the 128-feature Bayesian RFF model and a five-member MLP ensemble across four scenarios, three development seeds and random/PoH/CALDER rules (108 trajectories). PoH had more mean confirmed hits than CALDER in all twelve scenario/model combinations; the ensemble retained 153 fixed-budget convergence warnings, so these small checks do not establish model-family superiority or optimizer convergence. Separate [one-factor sensitivity runs](../results/reference_sensitivity/summary.md) contain 216 trajectories covering pilot fractions 0.1/0.2/0.3, batch sizes 5/10/20 and noise scales 0.05/0.2/0.5. Small positive CALDER differences occurred in some settings, for example C at pilot fraction 0.1 (+0.67 hits, exploratory CI [−2, 3]) and D at noise scale 0.5 (+1.0 [−2, 4]); these uncertain three-seed effects were retained. Pilot-size comparisons change both initial cost and the pilot-derived threshold, so they cannot be interpreted as isolated acquisition-rule effects across settings. No model or biological parameter was chosen from these outcomes. Full seed-level data and all group-coverage intervals are preserved in [ablations.csv](../results/reference_ablations/ablations.csv) and [paired_differences.csv](../results/reference_ablations/paired_differences.csv).

## 8. Limitations

Noisy pseudobulk response magnitudes depend on filtering, feature choice and cell sampling. Bootstrap uncertainty omits some experimental confounding. Shared controls can correlate estimation errors. Retrospective replay provides no prospective laboratory validation. Prequential examples are actively selected and therefore biased toward the acquisition distribution. Small pilot calibration reserves limit calibration power and can force fallback. Interval coverage is descriptive and has no exchangeability-based guarantee under adaptive selection.

Network annotations can be incomplete, inaccurate or context-insensitive. Feature clusters do not establish biological mechanisms. Conditional independence in coverage utility is approximate. RPE1 reruns evaluate a new context, whereas direct K562-to-RPE1 predictive transfer is a different optional experiment. Synthetic success cannot establish biological performance. A high response score does not establish therapeutic usefulness, and intervention effects do not automatically establish direct causal regulatory edges.

The harmonized scPerturb files already contain an upstream curated gene universe (8,563/8,749 variables). CALDER controls all newly fitted transformations, but cannot reconstruct every upstream expression-filtering choice from these files alone. Leakage prevention therefore applies to the delivered pipeline conditional on this released universe; it does not certify that the original release was prepared without access to intervention data. Original guide-level efficacy, cell-state and experimental-batch confounding are not fully corrected by the current pseudobulk score. The five artificial feature clusters are not independent biological-group annotations, limiting the biological interpretation of the central coverage question. Context-specific thresholds and different candidate pools also limit direct comparisons of absolute counts.

## 9. Discussion

The reference synthetic runs give a negative answer to the proposed joint-improvement hypothesis. Ordinary Probability-of-Hit discovered more measured and latent hits at the displayed central threshold. CALDER improved group coverage only in the separated-mode scenario B, at a measurable hit-count cost. In nearly saturated groups, a diversity term has little primary-outcome headroom. These outcomes reinforce the need to retain hit count and coverage as distinct endpoints.

Probability calibration also failed under this protocol: the pilot-reserved mapping worsened the displayed Brier/ECE diagnostics, and removing it improved C's discovery count. A fixed small reserve, finite-sample isotonic estimates and changing surrogate scores are plausible contributors. A monotone map can preserve single-candidate rankings while changing the scale and saturation of a batch diversity objective; better or worse probability metrics therefore need not translate mechanically into corresponding acquisition gains. These mechanisms require follow-up experiments rather than retrospective storytelling.

The completed biological replay strengthens the negative hit-discovery finding under the frozen reference method: PoH found more hits in both K562 essential and RPE1. K562 calibration improved despite poorer selection, while RPE1 showed little Brier improvement and marked interval undercoverage. Coverage findings concern artificial control-feature clusters and cannot establish wider pathway discovery. These are fixed retrospective pools with ten algorithmic seeds, not independent biological replications or prospective experimental validation.

## 10. Conclusion

CALDER provides a reproducible, leakage-tested framework with controlled synthetic experiments, component ablations, exploratory model/sensitivity checks and retrospective public-data K562/RPE1 validation under a frozen method. The default calibrated, diversified rule found fewer hits than ordinary Probability-of-Hit in the main synthetic landscapes and both biological contexts at the central threshold. Small proxy-group coverage gains came with reduced hit discovery, and better K562 probability metrics did not improve acquisition efficiency. These negative and mixed outcomes are the demonstrated findings. No real laboratory experiments were performed by CALDER; no therapeutic benefit, causal regulatory mechanism or universal policy advantage is established.

## 11. References

See the numbered [bibliography](bibliography.md).
