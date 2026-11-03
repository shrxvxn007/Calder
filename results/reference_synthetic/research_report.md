# CALDER: Uncertainty-Calibrated, Diversity-Aware Active Learning for Budgeted Genetic Perturbation Discovery

## Abstract

We evaluated budgeted perturbation selection using 1080 saved trajectories, 9 selection/reference rules and 10 random seeds in controlled synthetic simulation. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -8 to -1; group-coverage mean differences ranged from -0.1 to 0.7. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

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

Run: `reference_synthetic`. Input configuration and full provenance are saved in `config.json` and `experiment_manifest.json`.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean | auc_hits_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.7 | calder | 10 | 43 | 4.4 | 1358 |
| synthetic | A | 0.7 | greedy_mean | 10 | 50.5 | 3.7 | 1621 |
| synthetic | A | 0.7 | oracle_reference | 10 | 65.8 | 5 | 1836 |
| synthetic | A | 0.7 | poh_diversity | 10 | 47.9 | 4.2 | 1505 |
| synthetic | A | 0.7 | probability_of_hit | 10 | 50.5 | 3.7 | 1623 |
| synthetic | A | 0.7 | random | 10 | 33.6 | 3.8 | 1146 |
| synthetic | A | 0.7 | thompson | 10 | 47.4 | 3.9 | 1497 |
| synthetic | A | 0.7 | ucb | 10 | 50.1 | 3.7 | 1601 |
| synthetic | A | 0.7 | uncertainty | 10 | 33.2 | 4.2 | 1151 |
| synthetic | A | 0.8 | calder | 10 | 32.6 | 3.1 | 1030 |
| synthetic | A | 0.8 | greedy_mean | 10 | 40.5 | 3.1 | 1290 |
| synthetic | A | 0.8 | oracle_reference | 10 | 52.3 | 3.6 | 1580 |
| synthetic | A | 0.8 | poh_diversity | 10 | 38.1 | 3.3 | 1163 |
| synthetic | A | 0.8 | probability_of_hit | 10 | 40.6 | 3.1 | 1288 |
| synthetic | A | 0.8 | random | 10 | 23.2 | 3.2 | 809.4 |
| synthetic | A | 0.8 | thompson | 10 | 36.4 | 3.3 | 1148 |
| synthetic | A | 0.8 | ucb | 10 | 40.6 | 3.2 | 1272 |
| synthetic | A | 0.8 | uncertainty | 10 | 23.4 | 3.6 | 827 |
| synthetic | A | 0.9 | calder | 10 | 18.7 | 1.6 | 589.9 |
| synthetic | A | 0.9 | greedy_mean | 10 | 25.2 | 1.8 | 784.9 |
| synthetic | A | 0.9 | oracle_reference | 10 | 32 | 1.8 | 1127 |
| synthetic | A | 0.9 | poh_diversity | 10 | 22 | 1.7 | 638.4 |
| synthetic | A | 0.9 | probability_of_hit | 10 | 25.1 | 1.7 | 772.6 |
| synthetic | A | 0.9 | random | 10 | 12.8 | 1.3 | 442.1 |
| synthetic | A | 0.9 | thompson | 10 | 21.6 | 1.7 | 662.7 |
| synthetic | A | 0.9 | ucb | 10 | 24.7 | 1.7 | 755.2 |
| synthetic | A | 0.9 | uncertainty | 10 | 12.5 | 1.2 | 444.4 |
| synthetic | B | 0.7 | calder | 10 | 31.5 | 5 | 1098 |
| synthetic | B | 0.7 | greedy_mean | 10 | 35 | 4.9 | 1190 |
| synthetic | B | 0.7 | oracle_reference | 10 | 63.2 | 5 | 1800 |
| synthetic | B | 0.7 | poh_diversity | 10 | 34.7 | 5 | 1183 |
| synthetic | B | 0.7 | probability_of_hit | 10 | 35.2 | 4.8 | 1190 |
| synthetic | B | 0.7 | random | 10 | 31.4 | 4.9 | 1114 |
| synthetic | B | 0.7 | thompson | 10 | 33.4 | 4.9 | 1139 |
| synthetic | B | 0.7 | ucb | 10 | 34.6 | 4.9 | 1186 |
| synthetic | B | 0.7 | uncertainty | 10 | 28.8 | 5 | 1049 |
| synthetic | B | 0.8 | calder | 10 | 20.5 | 5 | 710.9 |
| synthetic | B | 0.8 | greedy_mean | 10 | 23.9 | 4.7 | 787.1 |
| synthetic | B | 0.8 | oracle_reference | 10 | 49.8 | 5 | 1463 |
| synthetic | B | 0.8 | poh_diversity | 10 | 24 | 4.8 | 787.9 |
| synthetic | B | 0.8 | probability_of_hit | 10 | 24.3 | 4.7 | 787.2 |
| synthetic | B | 0.8 | random | 10 | 20.3 | 4.8 | 720.9 |
| synthetic | B | 0.8 | thompson | 10 | 21.7 | 4.9 | 720 |
| synthetic | B | 0.8 | ucb | 10 | 23.1 | 4.7 | 762.3 |
| synthetic | B | 0.8 | uncertainty | 10 | 18.4 | 4.9 | 649 |
| synthetic | B | 0.9 | calder | 10 | 10.9 | 3.9 | 384.5 |
| synthetic | B | 0.9 | greedy_mean | 10 | 13.9 | 3.8 | 448.6 |
| synthetic | B | 0.9 | oracle_reference | 10 | 27.3 | 4.5 | 983.3 |
| synthetic | B | 0.9 | poh_diversity | 10 | 12.5 | 3.9 | 412.5 |
| synthetic | B | 0.9 | probability_of_hit | 10 | 14 | 3.8 | 456.8 |
| synthetic | B | 0.9 | random | 10 | 11.5 | 4.1 | 399.8 |
| synthetic | B | 0.9 | thompson | 10 | 11.9 | 3.9 | 394.8 |
| synthetic | B | 0.9 | ucb | 10 | 12.8 | 3.8 | 431.2 |
| synthetic | B | 0.9 | uncertainty | 10 | 10 | 3.8 | 348.9 |
| synthetic | C | 0.7 | calder | 10 | 33 | 5 | 1151 |
| synthetic | C | 0.7 | greedy_mean | 10 | 36.5 | 5 | 1237 |
| synthetic | C | 0.7 | oracle_reference | 10 | 64.6 | 5 | 1826 |
| synthetic | C | 0.7 | poh_diversity | 10 | 35.9 | 5 | 1221 |
| synthetic | C | 0.7 | probability_of_hit | 10 | 35.8 | 5 | 1224 |
| synthetic | C | 0.7 | random | 10 | 33.6 | 5 | 1176 |
| synthetic | C | 0.7 | thompson | 10 | 33.3 | 5 | 1162 |
| synthetic | C | 0.7 | ucb | 10 | 36 | 5 | 1219 |
| synthetic | C | 0.7 | uncertainty | 10 | 30.9 | 5 | 1100 |
| synthetic | C | 0.8 | calder | 10 | 23.6 | 4.9 | 828.5 |
| synthetic | C | 0.8 | greedy_mean | 10 | 27.7 | 5 | 919.6 |
| synthetic | C | 0.8 | oracle_reference | 10 | 57.6 | 5 | 1607 |
| synthetic | C | 0.8 | poh_diversity | 10 | 26.4 | 5 | 893.7 |
| synthetic | C | 0.8 | probability_of_hit | 10 | 26.7 | 5 | 906.5 |
| synthetic | C | 0.8 | random | 10 | 25.4 | 5 | 901.5 |
| synthetic | C | 0.8 | thompson | 10 | 24.8 | 5 | 860.7 |
| synthetic | C | 0.8 | ucb | 10 | 27.2 | 5 | 904.9 |
| synthetic | C | 0.8 | uncertainty | 10 | 23 | 5 | 813 |
| synthetic | C | 0.9 | calder | 10 | 16.1 | 4.9 | 578.2 |
| synthetic | C | 0.9 | greedy_mean | 10 | 19.4 | 4.6 | 634.1 |
| synthetic | C | 0.9 | oracle_reference | 10 | 43.4 | 5 | 1328 |
| synthetic | C | 0.9 | poh_diversity | 10 | 18.3 | 5 | 603.9 |
| synthetic | C | 0.9 | probability_of_hit | 10 | 18.9 | 4.6 | 631.6 |
| synthetic | C | 0.9 | random | 10 | 17.6 | 4.8 | 608.1 |
| synthetic | C | 0.9 | thompson | 10 | 16.7 | 4.8 | 579.7 |
| synthetic | C | 0.9 | ucb | 10 | 19 | 4.7 | 624.2 |
| synthetic | C | 0.9 | uncertainty | 10 | 15.8 | 4.7 | 546.6 |
| synthetic | D | 0.7 | calder | 10 | 32.1 | 5 | 1115 |
| synthetic | D | 0.7 | greedy_mean | 10 | 33.3 | 5 | 1140 |
| synthetic | D | 0.7 | oracle_reference | 10 | 64.1 | 5 | 1774 |
| synthetic | D | 0.7 | poh_diversity | 10 | 33.4 | 5 | 1142 |
| synthetic | D | 0.7 | probability_of_hit | 10 | 33.1 | 4.9 | 1153 |
| synthetic | D | 0.7 | random | 10 | 31.9 | 4.9 | 1089 |
| synthetic | D | 0.7 | thompson | 10 | 33.6 | 5 | 1142 |
| synthetic | D | 0.7 | ucb | 10 | 33 | 4.9 | 1118 |
| synthetic | D | 0.7 | uncertainty | 10 | 31.9 | 5 | 1102 |
| synthetic | D | 0.8 | calder | 10 | 22.5 | 5 | 802.6 |
| synthetic | D | 0.8 | greedy_mean | 10 | 23.9 | 5 | 833.8 |
| synthetic | D | 0.8 | oracle_reference | 10 | 52.2 | 5 | 1532 |
| synthetic | D | 0.8 | poh_diversity | 10 | 23.5 | 5 | 826 |
| synthetic | D | 0.8 | probability_of_hit | 10 | 23.8 | 5 | 841.1 |
| synthetic | D | 0.8 | random | 10 | 20.6 | 4.9 | 748.8 |
| synthetic | D | 0.8 | thompson | 10 | 24.4 | 4.8 | 842.3 |
| synthetic | D | 0.8 | ucb | 10 | 22.8 | 4.9 | 788.6 |
| synthetic | D | 0.8 | uncertainty | 10 | 22 | 4.8 | 769.5 |
| synthetic | D | 0.9 | calder | 10 | 13.1 | 4.9 | 486.4 |
| synthetic | D | 0.9 | greedy_mean | 10 | 15 | 4.7 | 527.6 |
| synthetic | D | 0.9 | oracle_reference | 10 | 31.2 | 5 | 1066 |
| synthetic | D | 0.9 | poh_diversity | 10 | 15.4 | 5 | 527.6 |
| synthetic | D | 0.9 | probability_of_hit | 10 | 15.1 | 4.8 | 529.7 |
| synthetic | D | 0.9 | random | 10 | 12.8 | 4.5 | 461.9 |
| synthetic | D | 0.9 | thompson | 10 | 14.4 | 4.6 | 496.4 |
| synthetic | D | 0.9 | ucb | 10 | 13.9 | 4.9 | 494.5 |
| synthetic | D | 0.9 | uncertainty | 10 | 13.1 | 4.6 | 459.8 |

Paired differences (policy minus reference), conditional on each fixed pool:

| dataset | scenario | threshold_quantile | policy | reference | metric | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.7 | calder | probability_of_hit | confirmed_hits | -7.5 | -9.7 | -5.5 | 10 |
| synthetic | A | 0.7 | calder | probability_of_hit | group_coverage | 0.7 | 0.2 | 1.2 | 10 |
| synthetic | A | 0.7 | calder | probability_of_hit | auc_hits | -265.7 | -328.4 | -209.4 | 10 |
| synthetic | A | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | 0 | -1 | 1.2 | 10 |
| synthetic | A | 0.7 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.7 | greedy_mean | probability_of_hit | auc_hits | -2.3 | -19.31 | 10.3 | 10 |
| synthetic | A | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | -2.6 | -4.2 | -1.1 | 10 |
| synthetic | A | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0.5 | 0.2 | 0.8 | 10 |
| synthetic | A | 0.7 | poh_diversity | probability_of_hit | auc_hits | -118.6 | -168.3 | -71.79 | 10 |
| synthetic | A | 0.7 | random | probability_of_hit | confirmed_hits | -16.9 | -18.4 | -15.3 | 10 |
| synthetic | A | 0.7 | random | probability_of_hit | group_coverage | 0.1 | -0.2 | 0.4 | 10 |
| synthetic | A | 0.7 | random | probability_of_hit | auc_hits | -477.1 | -535.6 | -419.2 | 10 |
| synthetic | A | 0.7 | thompson | probability_of_hit | confirmed_hits | -3.1 | -4 | -2.2 | 10 |
| synthetic | A | 0.7 | thompson | probability_of_hit | group_coverage | 0.2 | -0.2 | 0.6 | 10 |
| synthetic | A | 0.7 | thompson | probability_of_hit | auc_hits | -126.6 | -158.7 | -90.79 | 10 |
| synthetic | A | 0.7 | ucb | probability_of_hit | confirmed_hits | -0.4 | -1.4 | 0.5 | 10 |
| synthetic | A | 0.7 | ucb | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| synthetic | A | 0.7 | ucb | probability_of_hit | auc_hits | -22.2 | -42.1 | 1.305 | 10 |
| synthetic | A | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -17.3 | -19.2 | -15.2 | 10 |
| synthetic | A | 0.7 | uncertainty | probability_of_hit | group_coverage | 0.5 | 0 | 1.1 | 10 |
| synthetic | A | 0.7 | uncertainty | probability_of_hit | auc_hits | -472.6 | -528.1 | -406.4 | 10 |
| synthetic | A | 0.8 | calder | probability_of_hit | confirmed_hits | -8 | -10.3 | -5.798 | 10 |
| synthetic | A | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.8 | calder | probability_of_hit | auc_hits | -258 | -319.7 | -194.1 | 10 |
| synthetic | A | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | -0.1 | -0.6 | 0.4 | 10 |
| synthetic | A | 0.8 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.8 | greedy_mean | probability_of_hit | auc_hits | 1.8 | -4.6 | 8.8 | 10 |
| synthetic | A | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -2.5 | -3.8 | -1.3 | 10 |
| synthetic | A | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | A | 0.8 | poh_diversity | probability_of_hit | auc_hits | -125.2 | -158.1 | -90.2 | 10 |
| synthetic | A | 0.8 | random | probability_of_hit | confirmed_hits | -17.4 | -19.9 | -14.7 | 10 |
| synthetic | A | 0.8 | random | probability_of_hit | group_coverage | 0.1 | -0.2 | 0.4 | 10 |
| synthetic | A | 0.8 | random | probability_of_hit | auc_hits | -478.8 | -556.9 | -395.8 | 10 |
| synthetic | A | 0.8 | thompson | probability_of_hit | confirmed_hits | -4.2 | -5.8 | -2.4 | 10 |
| synthetic | A | 0.8 | thompson | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | A | 0.8 | thompson | probability_of_hit | auc_hits | -140.1 | -194.5 | -82.19 | 10 |
| synthetic | A | 0.8 | ucb | probability_of_hit | confirmed_hits | 0 | -0.7 | 0.8 | 10 |
| synthetic | A | 0.8 | ucb | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.8 | ucb | probability_of_hit | auc_hits | -16.6 | -38 | 9.9 | 10 |
| synthetic | A | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -17.2 | -19.1 | -15.1 | 10 |
| synthetic | A | 0.8 | uncertainty | probability_of_hit | group_coverage | 0.5 | 0.1 | 1 | 10 |
| synthetic | A | 0.8 | uncertainty | probability_of_hit | auc_hits | -461.2 | -532.8 | -383.5 | 10 |
| synthetic | A | 0.9 | calder | probability_of_hit | confirmed_hits | -6.4 | -8.7 | -4.1 | 10 |
| synthetic | A | 0.9 | calder | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | A | 0.9 | calder | probability_of_hit | auc_hits | -182.7 | -244.2 | -128.3 | 10 |
| synthetic | A | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | 0.1 | -0.8 | 1.4 | 10 |
| synthetic | A | 0.9 | greedy_mean | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.9 | greedy_mean | probability_of_hit | auc_hits | 12.3 | -11.8 | 40.21 | 10 |
| synthetic | A | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | -3.1 | -4.9 | -1.4 | 10 |
| synthetic | A | 0.9 | poh_diversity | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| synthetic | A | 0.9 | poh_diversity | probability_of_hit | auc_hits | -134.2 | -168.3 | -100.6 | 10 |
| synthetic | A | 0.9 | random | probability_of_hit | confirmed_hits | -12.3 | -14.5 | -9.9 | 10 |
| synthetic | A | 0.9 | random | probability_of_hit | group_coverage | -0.4 | -0.7 | -0.1 | 10 |
| synthetic | A | 0.9 | random | probability_of_hit | auc_hits | -330.5 | -383.4 | -265.6 | 10 |
| synthetic | A | 0.9 | thompson | probability_of_hit | confirmed_hits | -3.5 | -5.102 | -1.8 | 10 |
| synthetic | A | 0.9 | thompson | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| synthetic | A | 0.9 | thompson | probability_of_hit | auc_hits | -109.9 | -160 | -54.1 | 10 |
| synthetic | A | 0.9 | ucb | probability_of_hit | confirmed_hits | -0.4 | -1.3 | 0.4 | 10 |
| synthetic | A | 0.9 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.9 | ucb | probability_of_hit | auc_hits | -17.4 | -47.8 | 15.4 | 10 |
| synthetic | A | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -12.6 | -14.8 | -10.1 | 10 |
| synthetic | A | 0.9 | uncertainty | probability_of_hit | group_coverage | -0.5 | -0.8 | -0.2 | 10 |
| synthetic | A | 0.9 | uncertainty | probability_of_hit | auc_hits | -328.2 | -381.9 | -263.3 | 10 |
| synthetic | B | 0.7 | calder | probability_of_hit | confirmed_hits | -3.7 | -6.4 | -1.4 | 10 |
| synthetic | B | 0.7 | calder | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | B | 0.7 | calder | probability_of_hit | auc_hits | -92.1 | -148.6 | -37.7 | 10 |
| synthetic | B | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | -0.2 | -1.4 | 0.8 | 10 |
| synthetic | B | 0.7 | greedy_mean | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.7 | greedy_mean | probability_of_hit | auc_hits | -0.3 | -9.4 | 8.3 | 10 |
| synthetic | B | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | -0.5 | -1.5 | 0 | 10 |
| synthetic | B | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | B | 0.7 | poh_diversity | probability_of_hit | auc_hits | -6.9 | -17.11 | 0 | 10 |
| synthetic | B | 0.7 | random | probability_of_hit | confirmed_hits | -3.8 | -6.9 | -0.7 | 10 |
| synthetic | B | 0.7 | random | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.7 | random | probability_of_hit | auc_hits | -76.1 | -159.2 | 14.5 | 10 |
| synthetic | B | 0.7 | thompson | probability_of_hit | confirmed_hits | -1.8 | -4.7 | 0.9025 | 10 |
| synthetic | B | 0.7 | thompson | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.7 | thompson | probability_of_hit | auc_hits | -51.7 | -125.9 | 15.21 | 10 |
| synthetic | B | 0.7 | ucb | probability_of_hit | confirmed_hits | -0.6 | -3.1 | 1.8 | 10 |
| synthetic | B | 0.7 | ucb | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.7 | ucb | probability_of_hit | auc_hits | -4.1 | -41.8 | 40.8 | 10 |
| synthetic | B | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -6.4 | -9.7 | -3.298 | 10 |
| synthetic | B | 0.7 | uncertainty | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | B | 0.7 | uncertainty | probability_of_hit | auc_hits | -141.4 | -221.6 | -65.4 | 10 |
| synthetic | B | 0.8 | calder | probability_of_hit | confirmed_hits | -3.8 | -6.2 | -1.4 | 10 |
| synthetic | B | 0.8 | calder | probability_of_hit | group_coverage | 0.3 | 0.1 | 0.6 | 10 |
| synthetic | B | 0.8 | calder | probability_of_hit | auc_hits | -76.3 | -132.7 | -26.39 | 10 |
| synthetic | B | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | -0.4 | -1.1 | 0.2 | 10 |
| synthetic | B | 0.8 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.8 | greedy_mean | probability_of_hit | auc_hits | -0.1 | -24.4 | 17.4 | 10 |
| synthetic | B | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -0.3 | -1.1 | 0.5 | 10 |
| synthetic | B | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.8 | poh_diversity | probability_of_hit | auc_hits | 0.7 | -9.2 | 10.3 | 10 |
| synthetic | B | 0.8 | random | probability_of_hit | confirmed_hits | -4 | -6.5 | -1.7 | 10 |
| synthetic | B | 0.8 | random | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | B | 0.8 | random | probability_of_hit | auc_hits | -66.3 | -134.4 | -1.095 | 10 |
| synthetic | B | 0.8 | thompson | probability_of_hit | confirmed_hits | -2.6 | -4.7 | -0.3 | 10 |
| synthetic | B | 0.8 | thompson | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | B | 0.8 | thompson | probability_of_hit | auc_hits | -67.2 | -123.8 | -10.9 | 10 |
| synthetic | B | 0.8 | ucb | probability_of_hit | confirmed_hits | -1.2 | -3.1 | 0.7 | 10 |
| synthetic | B | 0.8 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.8 | ucb | probability_of_hit | auc_hits | -24.9 | -60.4 | 7.6 | 10 |
| synthetic | B | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -5.9 | -8.6 | -3.9 | 10 |
| synthetic | B | 0.8 | uncertainty | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | B | 0.8 | uncertainty | probability_of_hit | auc_hits | -138.2 | -213.1 | -74.29 | 10 |
| synthetic | B | 0.9 | calder | probability_of_hit | confirmed_hits | -3.1 | -4.4 | -1.9 | 10 |
| synthetic | B | 0.9 | calder | probability_of_hit | group_coverage | 0.1 | -0.2 | 0.4 | 10 |
| synthetic | B | 0.9 | calder | probability_of_hit | auc_hits | -72.3 | -100 | -43.4 | 10 |
| synthetic | B | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | -0.1 | -1 | 0.7025 | 10 |
| synthetic | B | 0.9 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.9 | greedy_mean | probability_of_hit | auc_hits | -8.2 | -27.1 | 5.402 | 10 |
| synthetic | B | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | -1.5 | -2.6 | -0.4 | 10 |
| synthetic | B | 0.9 | poh_diversity | probability_of_hit | group_coverage | 0.1 | -0.5 | 0.6 | 10 |
| synthetic | B | 0.9 | poh_diversity | probability_of_hit | auc_hits | -44.3 | -79.1 | -17.1 | 10 |
| synthetic | B | 0.9 | random | probability_of_hit | confirmed_hits | -2.5 | -4.3 | -0.8 | 10 |
| synthetic | B | 0.9 | random | probability_of_hit | group_coverage | 0.3 | -0.1 | 0.8 | 10 |
| synthetic | B | 0.9 | random | probability_of_hit | auc_hits | -57 | -103.1 | -11.7 | 10 |
| synthetic | B | 0.9 | thompson | probability_of_hit | confirmed_hits | -2.1 | -3.6 | -0.7 | 10 |
| synthetic | B | 0.9 | thompson | probability_of_hit | group_coverage | 0.1 | -0.4 | 0.6 | 10 |
| synthetic | B | 0.9 | thompson | probability_of_hit | auc_hits | -62 | -110.6 | -18.3 | 10 |
| synthetic | B | 0.9 | ucb | probability_of_hit | confirmed_hits | -1.2 | -2.4 | 0 | 10 |
| synthetic | B | 0.9 | ucb | probability_of_hit | group_coverage | 0 | -0.4 | 0.4 | 10 |
| synthetic | B | 0.9 | ucb | probability_of_hit | auc_hits | -25.6 | -50.3 | -1.3 | 10 |
| synthetic | B | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -4 | -5.602 | -2.3 | 10 |
| synthetic | B | 0.9 | uncertainty | probability_of_hit | group_coverage | 0 | -0.6 | 0.5025 | 10 |
| synthetic | B | 0.9 | uncertainty | probability_of_hit | auc_hits | -107.9 | -163.6 | -53 | 10 |
| synthetic | C | 0.7 | calder | probability_of_hit | confirmed_hits | -2.8 | -4.4 | -1 | 10 |
| synthetic | C | 0.7 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder | probability_of_hit | auc_hits | -73.6 | -115.9 | -32.6 | 10 |
| synthetic | C | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | 0.7 | -0.7 | 1.9 | 10 |
| synthetic | C | 0.7 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | greedy_mean | probability_of_hit | auc_hits | 12.4 | -2.5 | 32.2 | 10 |
| synthetic | C | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | poh_diversity | probability_of_hit | auc_hits | -3.4 | -10.2 | 0 | 10 |
| synthetic | C | 0.7 | random | probability_of_hit | confirmed_hits | -2.2 | -4.2 | -0.1 | 10 |
| synthetic | C | 0.7 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | random | probability_of_hit | auc_hits | -48.2 | -118.1 | 20.1 | 10 |
| synthetic | C | 0.7 | thompson | probability_of_hit | confirmed_hits | -2.5 | -4.1 | -1.1 | 10 |
| synthetic | C | 0.7 | thompson | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | thompson | probability_of_hit | auc_hits | -62.2 | -108.1 | -19.29 | 10 |
| synthetic | C | 0.7 | ucb | probability_of_hit | confirmed_hits | 0.2 | -1.002 | 1.6 | 10 |
| synthetic | C | 0.7 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | ucb | probability_of_hit | auc_hits | -5.5 | -34.2 | 19.6 | 10 |
| synthetic | C | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -4.9 | -6.8 | -2.9 | 10 |
| synthetic | C | 0.7 | uncertainty | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | uncertainty | probability_of_hit | auc_hits | -124.4 | -184.3 | -61 | 10 |
| synthetic | C | 0.8 | calder | probability_of_hit | confirmed_hits | -3.1 | -4.5 | -1.8 | 10 |
| synthetic | C | 0.8 | calder | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | C | 0.8 | calder | probability_of_hit | auc_hits | -78 | -106.4 | -48.3 | 10 |
| synthetic | C | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | 1 | -0.5 | 2.4 | 10 |
| synthetic | C | 0.8 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | greedy_mean | probability_of_hit | auc_hits | 13.1 | 0.3 | 28 | 10 |
| synthetic | C | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -0.3 | -0.9 | 0 | 10 |
| synthetic | C | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | poh_diversity | probability_of_hit | auc_hits | -12.8 | -31.6 | 0 | 10 |
| synthetic | C | 0.8 | random | probability_of_hit | confirmed_hits | -1.3 | -2.9 | 0.5 | 10 |
| synthetic | C | 0.8 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | random | probability_of_hit | auc_hits | -5 | -61.7 | 53 | 10 |
| synthetic | C | 0.8 | thompson | probability_of_hit | confirmed_hits | -1.9 | -3.6 | -0.5 | 10 |
| synthetic | C | 0.8 | thompson | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | thompson | probability_of_hit | auc_hits | -45.8 | -101 | 0.8075 | 10 |
| synthetic | C | 0.8 | ucb | probability_of_hit | confirmed_hits | 0.5 | -0.2 | 1.3 | 10 |
| synthetic | C | 0.8 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | ucb | probability_of_hit | auc_hits | -1.6 | -21.8 | 17.3 | 10 |
| synthetic | C | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -3.7 | -4.8 | -2.7 | 10 |
| synthetic | C | 0.8 | uncertainty | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | uncertainty | probability_of_hit | auc_hits | -93.5 | -133.5 | -52.5 | 10 |
| synthetic | C | 0.9 | calder | probability_of_hit | confirmed_hits | -2.8 | -3.802 | -1.8 | 10 |
| synthetic | C | 0.9 | calder | probability_of_hit | group_coverage | 0.3 | 0 | 0.6 | 10 |
| synthetic | C | 0.9 | calder | probability_of_hit | auc_hits | -53.4 | -83 | -23.3 | 10 |
| synthetic | C | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | 0.5 | -0.5 | 1.4 | 10 |
| synthetic | C | 0.9 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.9 | greedy_mean | probability_of_hit | auc_hits | 2.5 | -6.7 | 12.1 | 10 |
| synthetic | C | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | -0.6 | -1.7 | 0.3 | 10 |
| synthetic | C | 0.9 | poh_diversity | probability_of_hit | group_coverage | 0.4 | 0 | 0.8 | 10 |
| synthetic | C | 0.9 | poh_diversity | probability_of_hit | auc_hits | -27.7 | -55.9 | -3.598 | 10 |
| synthetic | C | 0.9 | random | probability_of_hit | confirmed_hits | -1.3 | -2.7 | -0.1 | 10 |
| synthetic | C | 0.9 | random | probability_of_hit | group_coverage | 0.2 | -0.2 | 0.7 | 10 |
| synthetic | C | 0.9 | random | probability_of_hit | auc_hits | -23.5 | -79.1 | 26.5 | 10 |
| synthetic | C | 0.9 | thompson | probability_of_hit | confirmed_hits | -2.2 | -3.8 | -0.4 | 10 |
| synthetic | C | 0.9 | thompson | probability_of_hit | group_coverage | 0.2 | -0.2 | 0.7 | 10 |
| synthetic | C | 0.9 | thompson | probability_of_hit | auc_hits | -51.9 | -110.1 | 3.5 | 10 |
| synthetic | C | 0.9 | ucb | probability_of_hit | confirmed_hits | 0.1 | -0.9 | 1.2 | 10 |
| synthetic | C | 0.9 | ucb | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.9 | ucb | probability_of_hit | auc_hits | -7.4 | -27.6 | 14 | 10 |
| synthetic | C | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -3.1 | -4 | -2.1 | 10 |
| synthetic | C | 0.9 | uncertainty | probability_of_hit | group_coverage | 0.1 | -0.2 | 0.4 | 10 |
| synthetic | C | 0.9 | uncertainty | probability_of_hit | auc_hits | -85 | -115.5 | -55.4 | 10 |
| synthetic | D | 0.7 | calder | probability_of_hit | confirmed_hits | -1 | -3.3 | 1.5 | 10 |
| synthetic | D | 0.7 | calder | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.7 | calder | probability_of_hit | auc_hits | -38 | -85.1 | 13 | 10 |
| synthetic | D | 0.7 | greedy_mean | probability_of_hit | confirmed_hits | 0.2 | -0.3 | 0.6 | 10 |
| synthetic | D | 0.7 | greedy_mean | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.7 | greedy_mean | probability_of_hit | auc_hits | -13.1 | -27 | -1.5 | 10 |
| synthetic | D | 0.7 | poh_diversity | probability_of_hit | confirmed_hits | 0.3 | 0 | 0.7 | 10 |
| synthetic | D | 0.7 | poh_diversity | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.7 | poh_diversity | probability_of_hit | auc_hits | -11.5 | -32.1 | 0 | 10 |
| synthetic | D | 0.7 | random | probability_of_hit | confirmed_hits | -1.2 | -2.8 | 0.4 | 10 |
| synthetic | D | 0.7 | random | probability_of_hit | group_coverage | 0 | -0.3 | 0.3 | 10 |
| synthetic | D | 0.7 | random | probability_of_hit | auc_hits | -64.5 | -133.6 | 3.325 | 10 |
| synthetic | D | 0.7 | thompson | probability_of_hit | confirmed_hits | 0.5 | -0.9 | 1.8 | 10 |
| synthetic | D | 0.7 | thompson | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.7 | thompson | probability_of_hit | auc_hits | -11 | -71.8 | 39.5 | 10 |
| synthetic | D | 0.7 | ucb | probability_of_hit | confirmed_hits | -0.1 | -1.2 | 0.9 | 10 |
| synthetic | D | 0.7 | ucb | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | ucb | probability_of_hit | auc_hits | -35.3 | -61.71 | -9.898 | 10 |
| synthetic | D | 0.7 | uncertainty | probability_of_hit | confirmed_hits | -1.2 | -3.5 | 1.1 | 10 |
| synthetic | D | 0.7 | uncertainty | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.7 | uncertainty | probability_of_hit | auc_hits | -51.6 | -124.9 | 21.02 | 10 |
| synthetic | D | 0.8 | calder | probability_of_hit | confirmed_hits | -1.3 | -3 | 0.5025 | 10 |
| synthetic | D | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | calder | probability_of_hit | auc_hits | -38.5 | -75.1 | -1.1 | 10 |
| synthetic | D | 0.8 | greedy_mean | probability_of_hit | confirmed_hits | 0.1 | -0.3 | 0.5 | 10 |
| synthetic | D | 0.8 | greedy_mean | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | greedy_mean | probability_of_hit | auc_hits | -7.3 | -17.8 | 1 | 10 |
| synthetic | D | 0.8 | poh_diversity | probability_of_hit | confirmed_hits | -0.3 | -0.9 | 0.2 | 10 |
| synthetic | D | 0.8 | poh_diversity | probability_of_hit | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | poh_diversity | probability_of_hit | auc_hits | -15.1 | -37.9 | 2.2 | 10 |
| synthetic | D | 0.8 | random | probability_of_hit | confirmed_hits | -3.2 | -5.202 | -1.098 | 10 |
| synthetic | D | 0.8 | random | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | D | 0.8 | random | probability_of_hit | auc_hits | -92.3 | -158.2 | -23.9 | 10 |
| synthetic | D | 0.8 | thompson | probability_of_hit | confirmed_hits | 0.6 | -1 | 2 | 10 |
| synthetic | D | 0.8 | thompson | probability_of_hit | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.8 | thompson | probability_of_hit | auc_hits | 1.2 | -54.1 | 46.8 | 10 |
| synthetic | D | 0.8 | ucb | probability_of_hit | confirmed_hits | -1 | -2.3 | 0.0025 | 10 |
| synthetic | D | 0.8 | ucb | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | D | 0.8 | ucb | probability_of_hit | auc_hits | -52.5 | -76.1 | -31.4 | 10 |
| synthetic | D | 0.8 | uncertainty | probability_of_hit | confirmed_hits | -1.8 | -4.3 | 0.9 | 10 |
| synthetic | D | 0.8 | uncertainty | probability_of_hit | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.8 | uncertainty | probability_of_hit | auc_hits | -71.6 | -136 | -4.185 | 10 |
| synthetic | D | 0.9 | calder | probability_of_hit | confirmed_hits | -2 | -3.5 | -0.6 | 10 |
| synthetic | D | 0.9 | calder | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.9 | calder | probability_of_hit | auc_hits | -43.3 | -87.6 | -1.393 | 10 |
| synthetic | D | 0.9 | greedy_mean | probability_of_hit | confirmed_hits | -0.1 | -1 | 0.7 | 10 |
| synthetic | D | 0.9 | greedy_mean | probability_of_hit | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | D | 0.9 | greedy_mean | probability_of_hit | auc_hits | -2.1 | -8.6 | 4.1 | 10 |
| synthetic | D | 0.9 | poh_diversity | probability_of_hit | confirmed_hits | 0.3 | 0.1 | 0.6 | 10 |
| synthetic | D | 0.9 | poh_diversity | probability_of_hit | group_coverage | 0.2 | 0 | 0.5 | 10 |
| synthetic | D | 0.9 | poh_diversity | probability_of_hit | auc_hits | -2.1 | -15.6 | 8.7 | 10 |
| synthetic | D | 0.9 | random | probability_of_hit | confirmed_hits | -2.3 | -4.1 | -0.4 | 10 |
| synthetic | D | 0.9 | random | probability_of_hit | group_coverage | -0.3 | -0.7 | 0.1 | 10 |
| synthetic | D | 0.9 | random | probability_of_hit | auc_hits | -67.8 | -127.1 | -12.4 | 10 |
| synthetic | D | 0.9 | thompson | probability_of_hit | confirmed_hits | -0.7 | -2.3 | 0.8 | 10 |
| synthetic | D | 0.9 | thompson | probability_of_hit | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.9 | thompson | probability_of_hit | auc_hits | -33.3 | -71.31 | 0.2 | 10 |
| synthetic | D | 0.9 | ucb | probability_of_hit | confirmed_hits | -1.2 | -2 | -0.4 | 10 |
| synthetic | D | 0.9 | ucb | probability_of_hit | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.9 | ucb | probability_of_hit | auc_hits | -35.2 | -60 | -14.6 | 10 |
| synthetic | D | 0.9 | uncertainty | probability_of_hit | confirmed_hits | -2 | -3.8 | -0.1 | 10 |
| synthetic | D | 0.9 | uncertainty | probability_of_hit | group_coverage | -0.2 | -0.5 | 0.2 | 10 |
| synthetic | D | 0.9 | uncertainty | probability_of_hit | auc_hits | -69.9 | -115.2 | -24.2 | 10 |

Full confidence intervals are in `results.csv`; seed endpoints are in `seed_level_results.csv`. Prequential latent/measured diagnostics remain separate in `calibration.csv`. Figure availability and links are recorded in [summary.md](summary.md).

## 7. Ablation studies

This run contains no completed ablation trajectories. No ablation effect is claimed. The predeclared ablation configuration must be executed separately.

## 8. Limitations

Noisy pseudobulk response magnitudes depend on filtering, feature choice and cell sampling. Bootstrap uncertainty omits some experimental confounding. Shared controls can correlate estimation errors. Retrospective replay provides no prospective laboratory validation. Prequential examples are actively selected and therefore biased toward the acquisition distribution. Small pilot calibration reserves limit calibration power and can force fallback. Interval coverage is descriptive and has no exchangeability-based guarantee under adaptive selection.

Network annotations can be incomplete, inaccurate or context-insensitive. Feature clusters do not establish biological mechanisms. Conditional independence in coverage utility is approximate. RPE1 reruns evaluate a new context, whereas direct K562-to-RPE1 predictive transfer is a different optional experiment. Synthetic success cannot establish biological performance. A high response score does not establish therapeutic usefulness, and intervention effects do not automatically establish direct causal regulatory edges.

## 9. Discussion

Interpret hit discovery and group coverage jointly without combining them into an arbitrary headline score. If diversity adds groups while reducing hits, quantify both changes. If calibration improves Brier score without improving acquisition, examine monotone ranking preservation, sample size, shifted score distributions and saturation of the coverage term. A more accurate model may still select fewer hits under a particular finite-budget rule.

## 10. Conclusion

We evaluated budgeted perturbation selection using 1080 saved trajectories, 9 selection/reference rules and 10 random seeds in controlled synthetic simulation. Across separately analyzed dataset/scenario/threshold strata, CALDER minus Probability-of-Hit confirmed-hit mean differences ranged from -8 to -1; group-coverage mean differences ranged from -0.1 to 0.7. These ranges describe separate effects rather than a pooled estimate. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

## 11. References

See the numbered [bibliography](bibliography.md).
