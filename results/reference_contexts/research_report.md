# CALDER: Uncertainty-Calibrated, Diversity-Aware Active Learning for Budgeted Genetic Perturbation Discovery

## Abstract

We evaluated budgeted perturbation selection using 540 saved trajectories, 9 selection/reference rules and 10 random seeds in retrospective biological replay. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -20.1 to -8.7; group-coverage mean differences ranged from 0 to 0.5. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority.

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

Measured-response predictive intervals add epistemic variance, estimated sampling variance and the configured model-discrepancy floor. The label `measured_predictive_model_assumed` identifies this likelihood assumption; it is heuristic for ridge/ensemble uncertainty. Latent-effect intervals use epistemic variance alone. Older saved runs recover the floor from their resolved configuration in memory while retaining original prediction CSVs.

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

## 6. Results

Run: `reference_contexts`. Input configuration and full provenance are saved in `config.json` and `experiment_manifest.json`.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean | auc_hits_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| replogle22k562 | replogle22k562 | 0.7 | calder | 10 | 146 | 5 | 3.118e+04 |
| replogle22k562 | replogle22k562 | 0.7 | greedy_mean | 10 | 157.3 | 4.3 | 3.381e+04 |
| replogle22k562 | replogle22k562 | 0.7 | oracle_reference | 10 | 399 | 5 | 7.738e+04 |
| replogle22k562 | replogle22k562 | 0.7 | poh_diversity | 10 | 162.5 | 4.7 | 3.434e+04 |
| replogle22k562 | replogle22k562 | 0.7 | probability_of_hit | 10 | 160 | 4.5 | 3.396e+04 |
| replogle22k562 | replogle22k562 | 0.7 | random | 10 | 136.9 | 4.5 | 2.977e+04 |
| replogle22k562 | replogle22k562 | 0.7 | thompson | 10 | 156.1 | 4.6 | 3.304e+04 |
| replogle22k562 | replogle22k562 | 0.7 | ucb | 10 | 161.8 | 4.4 | 3.397e+04 |
| replogle22k562 | replogle22k562 | 0.7 | uncertainty | 10 | 142.5 | 4.7 | 3.143e+04 |
| replogle22k562 | replogle22k562 | 0.8 | calder | 10 | 95.9 | 4.5 | 2.087e+04 |
| replogle22k562 | replogle22k562 | 0.8 | greedy_mean | 10 | 108.2 | 4.1 | 2.341e+04 |
| replogle22k562 | replogle22k562 | 0.8 | oracle_reference | 10 | 360.4 | 4.5 | 7.325e+04 |
| replogle22k562 | replogle22k562 | 0.8 | poh_diversity | 10 | 107.3 | 4.4 | 2.306e+04 |
| replogle22k562 | replogle22k562 | 0.8 | probability_of_hit | 10 | 110.2 | 4.2 | 2.346e+04 |
| replogle22k562 | replogle22k562 | 0.8 | random | 10 | 93.7 | 4.2 | 2.048e+04 |
| replogle22k562 | replogle22k562 | 0.8 | thompson | 10 | 107.8 | 4.2 | 2.266e+04 |
| replogle22k562 | replogle22k562 | 0.8 | ucb | 10 | 110.5 | 4.1 | 2.334e+04 |
| replogle22k562 | replogle22k562 | 0.8 | uncertainty | 10 | 96.4 | 4.5 | 2.127e+04 |
| replogle22k562 | replogle22k562 | 0.9 | calder | 10 | 54 | 4 | 1.12e+04 |
| replogle22k562 | replogle22k562 | 0.9 | greedy_mean | 10 | 62.6 | 4 | 1.327e+04 |
| replogle22k562 | replogle22k562 | 0.9 | oracle_reference | 10 | 231.7 | 4 | 5.862e+04 |
| replogle22k562 | replogle22k562 | 0.9 | poh_diversity | 10 | 66.6 | 4 | 1.36e+04 |
| replogle22k562 | replogle22k562 | 0.9 | probability_of_hit | 10 | 62.7 | 4 | 1.33e+04 |
| replogle22k562 | replogle22k562 | 0.9 | random | 10 | 58 | 4 | 1.249e+04 |
| replogle22k562 | replogle22k562 | 0.9 | thompson | 10 | 63.6 | 4 | 1.296e+04 |
| replogle22k562 | replogle22k562 | 0.9 | ucb | 10 | 64.4 | 4 | 1.335e+04 |
| replogle22k562 | replogle22k562 | 0.9 | uncertainty | 10 | 55.2 | 4 | 1.183e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | calder | 10 | 181.9 | 5 | 4.214e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | greedy_mean | 10 | 204.4 | 5 | 4.796e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | oracle_reference | 10 | 443.7 | 5 | 9.582e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | poh_diversity | 10 | 200.5 | 5 | 4.765e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | probability_of_hit | 10 | 200.7 | 5 | 4.772e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | random | 10 | 167.8 | 4.9 | 4.073e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | thompson | 10 | 200.9 | 4.9 | 4.692e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | ucb | 10 | 202.4 | 5 | 4.665e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | uncertainty | 10 | 165.4 | 5 | 4.065e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | calder | 10 | 126.5 | 5 | 2.877e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | greedy_mean | 10 | 147.9 | 4.8 | 3.44e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | oracle_reference | 10 | 401.8 | 5 | 8.976e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | poh_diversity | 10 | 144 | 4.8 | 3.414e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | probability_of_hit | 10 | 146.6 | 4.7 | 3.47e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | random | 10 | 117.6 | 4.6 | 2.87e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | thompson | 10 | 146.6 | 4.7 | 3.379e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | ucb | 10 | 147 | 4.8 | 3.316e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | uncertainty | 10 | 113.3 | 4.9 | 2.771e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | calder | 10 | 70.6 | 5 | 1.612e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | greedy_mean | 10 | 91.8 | 4.4 | 2.105e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | oracle_reference | 10 | 267.5 | 5 | 7.367e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | poh_diversity | 10 | 89.1 | 4.4 | 2.019e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | probability_of_hit | 10 | 89.2 | 4.7 | 2.054e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | random | 10 | 70.3 | 4.2 | 1.658e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | thompson | 10 | 90.7 | 4.5 | 2.064e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | ucb | 10 | 90.8 | 4.4 | 1.964e+04 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | uncertainty | 10 | 65.4 | 4.6 | 1.553e+04 |

Paired differences (policy minus reference), conditional on each fixed pool:

| dataset | scenario | threshold_quantile | policy | reference | metric | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| replogle22k562 | replogle22k562 | 0.7 | calder | probability_of_hit | confirmed_hits | -14 | -19.1 | -8.7 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | calder | probability_of_hit | group_coverage | 0.5 | 0.2 | 0.8 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | calder | probability_of_hit | auc_hits | -2781 | -3849 | -1641 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | -2.7 | -6.602 | 1.8 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | greedy_mean | probability_of_hit | group_coverage | -0.2 | -0.5 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | greedy_mean | probability_of_hit | auc_hits | -148.2 | -1129 | 995.5 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | 2.5 | -1.402 | 7.2 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | poh_diversity | probability_of_hit | auc_hits | 379.2 | -524.5 | 1638 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | random | probability_of_hit | confirmed_hits | -23.1 | -30.5 | -14.4 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | random | probability_of_hit | group_coverage | 0 | -0.4 | 0.4 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | random | probability_of_hit | auc_hits | -4190 | -6022 | -1918 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | thompson | probability_of_hit | confirmed_hits | -3.9 | -11 | 3.6 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | thompson | probability_of_hit | group_coverage | 0.1 | -0.2 | 0.4 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | thompson | probability_of_hit | auc_hits | -923.2 | -2863 | 994.4 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | ucb | probability_of_hit | confirmed_hits | 1.8 | -4.2 | 7.9 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | ucb | probability_of_hit | group_coverage | -0.1 | -0.4 | 0.2 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | ucb | probability_of_hit | auc_hits | 7.8 | -1612 | 1444 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -17.5 | -22.7 | -12.4 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | uncertainty | probability_of_hit | group_coverage | 0.2 | -0.2 | 0.6 | 10 |
| replogle22k562 | replogle22k562 | 0.7 | uncertainty | probability_of_hit | auc_hits | -2536 | -3549 | -1407 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | calder | probability_of_hit | confirmed_hits | -14.3 | -18 | -10.6 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | calder | probability_of_hit | group_coverage | 0.3 | 0.1 | 0.6 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | calder | probability_of_hit | auc_hits | -2595 | -3135 | -2041 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | -2 | -6.9 | 3.8 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | greedy_mean | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | greedy_mean | probability_of_hit | auc_hits | -55.8 | -720 | 511.2 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -2.9 | -6.6 | 0.9 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | poh_diversity | probability_of_hit | auc_hits | -408.5 | -854.9 | 25.53 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | random | probability_of_hit | confirmed_hits | -16.5 | -24.7 | -7.798 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | random | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | random | probability_of_hit | auc_hits | -2982 | -4475 | -1286 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | thompson | probability_of_hit | confirmed_hits | -2.4 | -6.6 | 2.202 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | thompson | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | thompson | probability_of_hit | auc_hits | -799.8 | -1531 | -40.94 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | ucb | probability_of_hit | confirmed_hits | 0.3 | -3.8 | 4.7 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | ucb | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | ucb | probability_of_hit | auc_hits | -123.6 | -681.8 | 501.7 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -13.8 | -19.2 | -8.998 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | uncertainty | probability_of_hit | group_coverage | 0.3 | 0.1 | 0.6 | 10 |
| replogle22k562 | replogle22k562 | 0.8 | uncertainty | probability_of_hit | auc_hits | -2197 | -2941 | -1437 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | calder | probability_of_hit | confirmed_hits | -8.7 | -13.7 | -3.8 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | calder | probability_of_hit | auc_hits | -2098 | -2994 | -1178 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | -0.1 | -3.5 | 3.8 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | greedy_mean | probability_of_hit | auc_hits | -33.6 | -608.7 | 556.5 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | 3.9 | 0.5 | 7.502 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | poh_diversity | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | poh_diversity | probability_of_hit | auc_hits | 298.1 | -144.6 | 689.3 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | random | probability_of_hit | confirmed_hits | -4.7 | -10.9 | 2.8 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | random | probability_of_hit | auc_hits | -817 | -1852 | 594.9 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | thompson | probability_of_hit | confirmed_hits | 0.9 | -3.1 | 4.802 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | thompson | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | thompson | probability_of_hit | auc_hits | -346.9 | -940.6 | 287 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | ucb | probability_of_hit | confirmed_hits | 1.7 | -1.2 | 4.602 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | ucb | probability_of_hit | auc_hits | 47 | -476.6 | 552.2 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -7.5 | -11.6 | -3.9 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | uncertainty | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22k562 | replogle22k562 | 0.9 | uncertainty | probability_of_hit | auc_hits | -1468 | -2211 | -768.9 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | calder | probability_of_hit | confirmed_hits | -18.8 | -24.1 | -13.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | calder | probability_of_hit | auc_hits | -5582 | -7277 | -3595 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | 3.7 | -2.2 | 9.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | greedy_mean | probability_of_hit | auc_hits | 239.5 | -806 | 1326 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | -0.2 | -4.802 | 4.4 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | poh_diversity | probability_of_hit | auc_hits | -72.3 | -878.3 | 842.8 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | random | probability_of_hit | confirmed_hits | -32.9 | -37.9 | -28.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | random | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | random | probability_of_hit | auc_hits | -6992 | -8005 | -5884 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | thompson | probability_of_hit | confirmed_hits | 0.2 | -3.6 | 4 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | thompson | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | thompson | probability_of_hit | auc_hits | -798 | -2090 | 632.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | ucb | probability_of_hit | confirmed_hits | 1.7 | -2.902 | 6.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | ucb | probability_of_hit | auc_hits | -1064 | -1898 | -242 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -35.3 | -39.4 | -31.1 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | uncertainty | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.7 | uncertainty | probability_of_hit | auc_hits | -7066 | -8405 | -5572 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | calder | probability_of_hit | confirmed_hits | -20.1 | -28 | -10.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | calder | probability_of_hit | group_coverage | 0.3 | 0 | 0.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | calder | probability_of_hit | auc_hits | -5925 | -7202 | -4243 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | 1.3 | -3.2 | 5.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | greedy_mean | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | greedy_mean | probability_of_hit | auc_hits | -294.8 | -1109 | 440.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -2.6 | -6.4 | 1.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | poh_diversity | probability_of_hit | auc_hits | -561.2 | -1879 | 719.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | random | probability_of_hit | confirmed_hits | -29 | -36.3 | -21.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | random | probability_of_hit | group_coverage | -0.1 | -0.5 | 0.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | random | probability_of_hit | auc_hits | -6000 | -7640 | -3950 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | thompson | probability_of_hit | confirmed_hits | 0 | -3.4 | 3.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | thompson | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | thompson | probability_of_hit | auc_hits | -905.1 | -2127 | 371.1 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | ucb | probability_of_hit | confirmed_hits | 0.4 | -3.9 | 4.4 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | ucb | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | ucb | probability_of_hit | auc_hits | -1534 | -2788 | -361.3 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -33.3 | -38.3 | -27.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | uncertainty | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.8 | uncertainty | probability_of_hit | auc_hits | -6986 | -8220 | -5509 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | calder | probability_of_hit | confirmed_hits | -18.6 | -26 | -11.1 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | calder | probability_of_hit | group_coverage | 0.3 | 0 | 0.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | calder | probability_of_hit | auc_hits | -4424 | -5939 | -2799 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | 2.6 | -0.7 | 5.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | greedy_mean | probability_of_hit | group_coverage | -0.3 | -0.6 | -0.1 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | greedy_mean | probability_of_hit | auc_hits | 503.7 | -44.14 | 1116 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | -0.1 | -3.6 | 2.9 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | poh_diversity | probability_of_hit | group_coverage | -0.3 | -0.6 | -0.1 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | poh_diversity | probability_of_hit | auc_hits | -357.8 | -841.5 | 121.7 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | random | probability_of_hit | confirmed_hits | -18.9 | -23.5 | -13.6 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | random | probability_of_hit | group_coverage | -0.5 | -0.8 | -0.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | random | probability_of_hit | auc_hits | -3966 | -5431 | -2350 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | thompson | probability_of_hit | confirmed_hits | 1.5 | -1.6 | 4.102 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | thompson | probability_of_hit | group_coverage | -0.2 | -0.6 | 0.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | thompson | probability_of_hit | auc_hits | 98.4 | -759.9 | 1017 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | ucb | probability_of_hit | confirmed_hits | 1.6 | -0.5 | 3.8 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | ucb | probability_of_hit | group_coverage | -0.3 | -0.6 | 0 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | ucb | probability_of_hit | auc_hits | -901.6 | -1584 | -247.5 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -23.8 | -28.7 | -18.8 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | uncertainty | probability_of_hit | group_coverage | -0.1 | -0.4 | 0.2 | 10 |
| replogle22rpe1 | replogle22rpe1 | 0.9 | uncertainty | probability_of_hit | auc_hits | -5015 | -6100 | -3930 | 10 |

Full confidence intervals are in `results.csv`; seed endpoints are in `seed_level_results.csv`. Prequential latent/measured diagnostics remain separate in `calibration.csv`. Figure availability and links are recorded in [summary.md](summary.md).

## 7. Ablation studies

This run contains no completed ablation trajectories. No ablation effect is claimed. The predeclared ablation configuration must be executed separately.

## 8. Limitations

Noisy pseudobulk response magnitudes depend on filtering, feature choice and cell sampling. Bootstrap uncertainty omits some experimental confounding. Shared controls can correlate estimation errors. Retrospective replay provides no prospective laboratory validation. Prequential examples are actively selected and therefore biased toward the acquisition distribution. Small pilot calibration reserves limit calibration power and can force fallback. Interval coverage is descriptive and has no exchangeability-based guarantee under adaptive selection.

Network annotations can be incomplete, inaccurate or context-insensitive. Feature clusters do not establish biological mechanisms. Conditional independence in coverage utility is approximate. RPE1 reruns evaluate a new context, whereas direct K562-to-RPE1 predictive transfer is a different optional experiment. Synthetic success cannot establish biological performance. A high response score does not establish therapeutic usefulness, and intervention effects do not automatically establish direct causal regulatory edges.

## 9. Discussion

Interpret hit discovery and group coverage jointly without combining them into an arbitrary headline score. If diversity adds groups while reducing hits, quantify both changes. If calibration improves Brier score without improving acquisition, examine monotone ranking preservation, sample size, shifted score distributions and saturation of the coverage term. A more accurate model may still select fewer hits under a particular finite-budget rule.

## 10. Conclusion

We evaluated budgeted perturbation selection using 540 saved trajectories, 9 selection/reference rules and 10 random seeds in retrospective biological replay. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -20.1 to -8.7; group-coverage mean differences ranged from 0 to 0.5. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority.

## 11. References

See the numbered [bibliography](bibliography.md).
