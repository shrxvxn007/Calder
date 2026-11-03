# CALDER: Uncertainty-Calibrated, Diversity-Aware Active Learning for Budgeted Genetic Perturbation Discovery

## Abstract

We evaluated budgeted perturbation selection using 36 saved trajectories, 3 selection/reference rules and 3 random seeds in controlled synthetic simulation. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -6 to -2.67; group-coverage mean differences ranged from 0 to 0.667. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

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

The primary configuration hash is `a76b99c6129af0d9adc3ed4cb1beaf456ba1beb0fa75f28503d87191a1ec2a52`; the ablation configuration hash is `ffe64235d97614f80d76e8a075a8d1ea29c7bccd5d185ddf3747e27cba281a81`. Both manifests record Git `15de4bebfb212eca645a74806bafbee35beb3a5b` with `git_dirty=true`; the commit alone therefore does not fully identify the execution working tree. Retain the delivered source alongside resolved configurations and manifests. The primary `seed_level_results.csv` SHA-256 is `8012c6cb4afda2aa45a27c64dfd965179621b7665e85233fe2d486671d0651f8`, and the ablation endpoint-table SHA-256 is `6a43bec664085c5b587ee5ee018c26814bbc696e9d73db3bfda4bdca86aa2821`. These four key saved files (configuration, seed endpoints, paired differences and calibration diagnostics) were checked against each run manifest. There were no repeated candidate indices within any saved prequential trajectory.

Primary acquisition/evaluation runtime was 29.67 seconds; the ablation run was 27.14 seconds on the recorded CPU environment, excluding subsequent figure/report generation. Process lifetime resident high-water marks were 492.3 MB and 443.9 MB respectively; these include imports and are not incremental per-trajectory allocations.

## 6. Results

Run: `ensemble`. Input configuration and full provenance are saved in `config.json` and `experiment_manifest.json`.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean | auc_hits_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.8 | calder | 3 | 42.67 | 2.667 | 1329 |
| synthetic | A | 0.8 | probability_of_hit | 3 | 46 | 2.667 | 1475 |
| synthetic | A | 0.8 | random | 3 | 22.67 | 3 | 810.3 |
| synthetic | B | 0.8 | calder | 3 | 20 | 5 | 716 |
| synthetic | B | 0.8 | probability_of_hit | 3 | 26 | 4.333 | 853 |
| synthetic | B | 0.8 | random | 3 | 20 | 5 | 740.3 |
| synthetic | C | 0.8 | calder | 3 | 21.67 | 5 | 734.7 |
| synthetic | C | 0.8 | probability_of_hit | 3 | 26 | 5 | 837.7 |
| synthetic | C | 0.8 | random | 3 | 24.67 | 5 | 860.3 |
| synthetic | D | 0.8 | calder | 3 | 23.67 | 5 | 780.7 |
| synthetic | D | 0.8 | probability_of_hit | 3 | 26.33 | 4.667 | 842.3 |
| synthetic | D | 0.8 | random | 3 | 19 | 5 | 694.3 |

Paired differences (policy minus reference), conditional on each fixed pool:

| dataset | scenario | threshold_quantile | policy | reference | metric | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.8 | calder | probability_of_hit | confirmed_hits | -3.333 | -8 | 0 | 3 |
| synthetic | A | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | A | 0.8 | calder | probability_of_hit | auc_hits | -146 | -347 | -31 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | confirmed_hits | -23.33 | -27 | -20 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | auc_hits | -665 | -682 | -650 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | confirmed_hits | -6 | -11 | -3 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | group_coverage | 0.6667 | 0 | 2 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | auc_hits | -137 | -201 | -80 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | confirmed_hits | -6 | -8 | -2 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | group_coverage | 0.6667 | 0 | 2 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | auc_hits | -112.7 | -171 | -7 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | confirmed_hits | -4.333 | -8 | -2 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | auc_hits | -103 | -191 | -53 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | confirmed_hits | -1.333 | -4 | 1 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | auc_hits | 22.67 | -84 | 105 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | confirmed_hits | -2.667 | -6 | -1 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | auc_hits | -61.67 | -138 | 37 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | confirmed_hits | -7.333 | -9 | -5 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | auc_hits | -148 | -166 | -116 | 3 |

Full confidence intervals are in `results.csv`; seed endpoints are in `seed_level_results.csv`. Prequential latent/measured diagnostics remain separate in `calibration.csv`. Figure availability and links are recorded in [summary.md](summary.md).

## 7. Ablation studies

This run contains no completed ablation trajectories. No ablation effect is claimed. The predeclared ablation configuration must be executed separately.

## 8. Limitations

Noisy pseudobulk response magnitudes depend on filtering, feature choice and cell sampling. Bootstrap uncertainty omits some experimental confounding. Shared controls can correlate estimation errors. Retrospective replay provides no prospective laboratory validation. Prequential examples are actively selected and therefore biased toward the acquisition distribution. Small pilot calibration reserves limit calibration power and can force fallback. Interval coverage is descriptive and has no exchangeability-based guarantee under adaptive selection.

Network annotations can be incomplete, inaccurate or context-insensitive. Feature clusters do not establish biological mechanisms. Conditional independence in coverage utility is approximate. RPE1 reruns evaluate a new context, whereas direct K562-to-RPE1 predictive transfer is a different optional experiment. Synthetic success cannot establish biological performance. A high response score does not establish therapeutic usefulness, and intervention effects do not automatically establish direct causal regulatory edges.

## 9. Discussion

The reference synthetic runs give a negative answer to the proposed joint-improvement hypothesis. Ordinary Probability-of-Hit discovered more measured and latent hits at the displayed central threshold. CALDER improved group coverage only in the separated-mode scenario B, at a measurable hit-count cost. In nearly saturated groups, a diversity term has little primary-outcome headroom. These outcomes reinforce the need to retain hit count and coverage as distinct endpoints.

Probability calibration also failed under this protocol: the pilot-reserved mapping worsened the displayed Brier/ECE diagnostics, and removing it improved C's discovery count. A fixed small reserve, finite-sample isotonic estimates and changing surrogate scores are plausible contributors. A monotone map can preserve single-candidate rankings while changing the scale and saturation of a batch diversity objective; better or worse probability metrics therefore need not translate mechanically into corresponding acquisition gains. These mechanisms require follow-up experiments rather than retrospective storytelling.

The findings concern four particular seeded synthetic pools, not population-wide claims about genetic perturbations. Ten seeds share each underlying pool. Biological evaluation must be performed with frozen methods, explicit source provenance and control-only representations before extending the interpretation to K562 or RPE1.

## 10. Conclusion

We evaluated budgeted perturbation selection using 36 saved trajectories, 3 selection/reference rules and 3 random seeds in controlled synthetic simulation. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -6 to -2.67; group-coverage mean differences ranged from 0 to 0.667. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

## 11. References

See the numbered [bibliography](bibliography.md).
