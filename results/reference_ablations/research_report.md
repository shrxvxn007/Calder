# CALDER: Uncertainty-Calibrated, Diversity-Aware Active Learning for Budgeted Genetic Perturbation Discovery

## Abstract

We evaluated budgeted perturbation selection using 840 saved trajectories, 7 selection/reference rules and 10 random seeds in controlled synthetic simulation. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

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

Run: `reference_ablations`. Input configuration and full provenance are saved in `config.json` and `experiment_manifest.json`.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean | auc_hits_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.7 | calder | 10 | 43 | 4.4 | 1358 |
| synthetic | A | 0.7 | calder_no_calibration | 10 | 43.4 | 4.2 | 1342 |
| synthetic | A | 0.7 | calder_no_diversity | 10 | 44.4 | 3.7 | 1389 |
| synthetic | A | 0.7 | calder_no_similarity | 10 | 45.9 | 4.1 | 1468 |
| synthetic | A | 0.7 | calder_no_uncertainty | 10 | 41.6 | 4.1 | 1348 |
| synthetic | A | 0.7 | calder_random_groups | 10 | 44.4 | 3.7 | 1389 |
| synthetic | A | 0.7 | calder_uninformative_features | 10 | 29.5 | 4.6 | 1080 |
| synthetic | A | 0.8 | calder | 10 | 32.6 | 3.1 | 1030 |
| synthetic | A | 0.8 | calder_no_calibration | 10 | 32.1 | 3.2 | 988.1 |
| synthetic | A | 0.8 | calder_no_diversity | 10 | 33.2 | 3.2 | 1055 |
| synthetic | A | 0.8 | calder_no_similarity | 10 | 36.3 | 3.2 | 1139 |
| synthetic | A | 0.8 | calder_no_uncertainty | 10 | 32.5 | 3.1 | 1033 |
| synthetic | A | 0.8 | calder_random_groups | 10 | 32.8 | 3.2 | 1054 |
| synthetic | A | 0.8 | calder_uninformative_features | 10 | 17.7 | 3.2 | 722.4 |
| synthetic | A | 0.9 | calder | 10 | 18.7 | 1.6 | 589.9 |
| synthetic | A | 0.9 | calder_no_calibration | 10 | 16.9 | 1.6 | 525.2 |
| synthetic | A | 0.9 | calder_no_diversity | 10 | 19.6 | 1.6 | 635.9 |
| synthetic | A | 0.9 | calder_no_similarity | 10 | 20.2 | 1.7 | 662.6 |
| synthetic | A | 0.9 | calder_no_uncertainty | 10 | 18.6 | 1.7 | 620.1 |
| synthetic | A | 0.9 | calder_random_groups | 10 | 20.2 | 1.6 | 653.4 |
| synthetic | A | 0.9 | calder_uninformative_features | 10 | 9.5 | 1.4 | 392.4 |
| synthetic | B | 0.7 | calder | 10 | 31.5 | 5 | 1098 |
| synthetic | B | 0.7 | calder_no_calibration | 10 | 32.9 | 5 | 1117 |
| synthetic | B | 0.7 | calder_no_diversity | 10 | 31.5 | 4.8 | 1099 |
| synthetic | B | 0.7 | calder_no_similarity | 10 | 32.4 | 5 | 1152 |
| synthetic | B | 0.7 | calder_no_uncertainty | 10 | 31.3 | 5 | 1095 |
| synthetic | B | 0.7 | calder_random_groups | 10 | 31.4 | 4.8 | 1104 |
| synthetic | B | 0.7 | calder_uninformative_features | 10 | 32.8 | 5 | 1142 |
| synthetic | B | 0.8 | calder | 10 | 20.5 | 5 | 710.9 |
| synthetic | B | 0.8 | calder_no_calibration | 10 | 21.6 | 4.9 | 728.5 |
| synthetic | B | 0.8 | calder_no_diversity | 10 | 20.2 | 4.6 | 703.8 |
| synthetic | B | 0.8 | calder_no_similarity | 10 | 21.8 | 4.7 | 731.5 |
| synthetic | B | 0.8 | calder_no_uncertainty | 10 | 20.7 | 4.7 | 716.4 |
| synthetic | B | 0.8 | calder_random_groups | 10 | 20.7 | 4.7 | 703.3 |
| synthetic | B | 0.8 | calder_uninformative_features | 10 | 19.1 | 4.9 | 671 |
| synthetic | B | 0.9 | calder | 10 | 10.9 | 3.9 | 384.5 |
| synthetic | B | 0.9 | calder_no_calibration | 10 | 11.1 | 3.9 | 384.7 |
| synthetic | B | 0.9 | calder_no_diversity | 10 | 11.4 | 3.7 | 407.1 |
| synthetic | B | 0.9 | calder_no_similarity | 10 | 12.2 | 4.2 | 395.5 |
| synthetic | B | 0.9 | calder_no_uncertainty | 10 | 10.7 | 3.9 | 386.9 |
| synthetic | B | 0.9 | calder_random_groups | 10 | 10.6 | 3.8 | 381 |
| synthetic | B | 0.9 | calder_uninformative_features | 10 | 9.5 | 4 | 342 |
| synthetic | C | 0.7 | calder | 10 | 33 | 5 | 1151 |
| synthetic | C | 0.7 | calder_no_calibration | 10 | 34.4 | 5 | 1186 |
| synthetic | C | 0.7 | calder_no_diversity | 10 | 33.1 | 5 | 1148 |
| synthetic | C | 0.7 | calder_no_similarity | 10 | 33.1 | 5 | 1158 |
| synthetic | C | 0.7 | calder_no_uncertainty | 10 | 33.7 | 5 | 1176 |
| synthetic | C | 0.7 | calder_random_groups | 10 | 33.2 | 5 | 1156 |
| synthetic | C | 0.7 | calder_uninformative_features | 10 | 31.7 | 5 | 1109 |
| synthetic | C | 0.8 | calder | 10 | 23.6 | 4.9 | 828.5 |
| synthetic | C | 0.8 | calder_no_calibration | 10 | 25.7 | 5 | 859.6 |
| synthetic | C | 0.8 | calder_no_diversity | 10 | 24 | 4.9 | 831.1 |
| synthetic | C | 0.8 | calder_no_similarity | 10 | 24.5 | 5 | 864.9 |
| synthetic | C | 0.8 | calder_no_uncertainty | 10 | 24.6 | 4.9 | 870.5 |
| synthetic | C | 0.8 | calder_random_groups | 10 | 23.7 | 4.9 | 830.3 |
| synthetic | C | 0.8 | calder_uninformative_features | 10 | 24 | 5 | 823.6 |
| synthetic | C | 0.9 | calder | 10 | 16.1 | 4.9 | 578.2 |
| synthetic | C | 0.9 | calder_no_calibration | 10 | 17 | 5 | 582.4 |
| synthetic | C | 0.9 | calder_no_diversity | 10 | 16.5 | 4.8 | 579.7 |
| synthetic | C | 0.9 | calder_no_similarity | 10 | 17.2 | 5 | 604.2 |
| synthetic | C | 0.9 | calder_no_uncertainty | 10 | 16.6 | 4.9 | 592.4 |
| synthetic | C | 0.9 | calder_random_groups | 10 | 16.6 | 4.7 | 593.3 |
| synthetic | C | 0.9 | calder_uninformative_features | 10 | 15.7 | 5 | 546.3 |
| synthetic | D | 0.7 | calder | 10 | 32.1 | 5 | 1115 |
| synthetic | D | 0.7 | calder_no_calibration | 10 | 31.8 | 5 | 1114 |
| synthetic | D | 0.7 | calder_no_diversity | 10 | 31.9 | 5 | 1109 |
| synthetic | D | 0.7 | calder_no_similarity | 10 | 32.6 | 5 | 1110 |
| synthetic | D | 0.7 | calder_no_uncertainty | 10 | 32 | 5 | 1109 |
| synthetic | D | 0.7 | calder_random_groups | 10 | 31.9 | 5 | 1109 |
| synthetic | D | 0.7 | calder_uninformative_features | 10 | 31.7 | 5 | 1086 |
| synthetic | D | 0.8 | calder | 10 | 22.5 | 5 | 802.6 |
| synthetic | D | 0.8 | calder_no_calibration | 10 | 23.2 | 5 | 819 |
| synthetic | D | 0.8 | calder_no_diversity | 10 | 21.7 | 4.9 | 781.4 |
| synthetic | D | 0.8 | calder_no_similarity | 10 | 22.8 | 5 | 805.5 |
| synthetic | D | 0.8 | calder_no_uncertainty | 10 | 22.9 | 5 | 814.4 |
| synthetic | D | 0.8 | calder_random_groups | 10 | 21.8 | 4.8 | 788.1 |
| synthetic | D | 0.8 | calder_uninformative_features | 10 | 21.4 | 5 | 760.7 |
| synthetic | D | 0.9 | calder | 10 | 13.1 | 4.9 | 486.4 |
| synthetic | D | 0.9 | calder_no_calibration | 10 | 13.9 | 5 | 495.9 |
| synthetic | D | 0.9 | calder_no_diversity | 10 | 12.8 | 4.7 | 460.3 |
| synthetic | D | 0.9 | calder_no_similarity | 10 | 13.9 | 4.9 | 486.7 |
| synthetic | D | 0.9 | calder_no_uncertainty | 10 | 13.1 | 4.9 | 490.2 |
| synthetic | D | 0.9 | calder_random_groups | 10 | 12.2 | 4.7 | 464.4 |
| synthetic | D | 0.9 | calder_uninformative_features | 10 | 12.5 | 4.7 | 460.2 |

Paired differences (policy minus reference), conditional on each fixed pool:

| dataset | scenario | threshold_quantile | policy | reference | metric | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.7 | calder_no_calibration | calder | confirmed_hits | 0.4 | -1.5 | 2.3 | 10 |
| synthetic | A | 0.7 | calder_no_calibration | calder | group_coverage | -0.2 | -0.8 | 0.3 | 10 |
| synthetic | A | 0.7 | calder_no_calibration | calder | auc_hits | -15.7 | -65.41 | 32.22 | 10 |
| synthetic | A | 0.7 | calder_no_diversity | calder | confirmed_hits | 1.4 | 0.3 | 2.4 | 10 |
| synthetic | A | 0.7 | calder_no_diversity | calder | group_coverage | -0.7 | -1.2 | -0.2 | 10 |
| synthetic | A | 0.7 | calder_no_diversity | calder | auc_hits | 31.8 | 13.5 | 48.1 | 10 |
| synthetic | A | 0.7 | calder_no_similarity | calder | confirmed_hits | 2.9 | 1.4 | 4.5 | 10 |
| synthetic | A | 0.7 | calder_no_similarity | calder | group_coverage | -0.3 | -0.8 | 0 | 10 |
| synthetic | A | 0.7 | calder_no_similarity | calder | auc_hits | 111 | 73.59 | 149.7 | 10 |
| synthetic | A | 0.7 | calder_no_uncertainty | calder | confirmed_hits | -1.4 | -3.2 | 0.2 | 10 |
| synthetic | A | 0.7 | calder_no_uncertainty | calder | group_coverage | -0.3 | -0.8 | 0 | 10 |
| synthetic | A | 0.7 | calder_no_uncertainty | calder | auc_hits | -9.2 | -50.2 | 29.5 | 10 |
| synthetic | A | 0.7 | calder_random_groups | calder | confirmed_hits | 1.4 | 0.3 | 2.4 | 10 |
| synthetic | A | 0.7 | calder_random_groups | calder | group_coverage | -0.7 | -1.2 | -0.2 | 10 |
| synthetic | A | 0.7 | calder_random_groups | calder | auc_hits | 31.8 | 13.5 | 48.1 | 10 |
| synthetic | A | 0.7 | calder_uninformative_features | calder | confirmed_hits | -13.5 | -16.6 | -10.4 | 10 |
| synthetic | A | 0.7 | calder_uninformative_features | calder | group_coverage | 0.2 | -0.3 | 0.7025 | 10 |
| synthetic | A | 0.7 | calder_uninformative_features | calder | auc_hits | -277.4 | -388 | -166.4 | 10 |
| synthetic | A | 0.8 | calder_no_calibration | calder | confirmed_hits | -0.5 | -2 | 1.2 | 10 |
| synthetic | A | 0.8 | calder_no_calibration | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.8 | calder_no_calibration | calder | auc_hits | -42.1 | -94.31 | 7.802 | 10 |
| synthetic | A | 0.8 | calder_no_diversity | calder | confirmed_hits | 0.6 | -0.5 | 1.9 | 10 |
| synthetic | A | 0.8 | calder_no_diversity | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.8 | calder_no_diversity | calder | auc_hits | 25.1 | -7.207 | 62.11 | 10 |
| synthetic | A | 0.8 | calder_no_similarity | calder | confirmed_hits | 3.7 | 0.9 | 6.3 | 10 |
| synthetic | A | 0.8 | calder_no_similarity | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.8 | calder_no_similarity | calder | auc_hits | 108.7 | 57.6 | 157.3 | 10 |
| synthetic | A | 0.8 | calder_no_uncertainty | calder | confirmed_hits | -0.1 | -3 | 2.8 | 10 |
| synthetic | A | 0.8 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.8 | calder_no_uncertainty | calder | auc_hits | 3 | -40.1 | 49.2 | 10 |
| synthetic | A | 0.8 | calder_random_groups | calder | confirmed_hits | 0.2 | -1.2 | 1.8 | 10 |
| synthetic | A | 0.8 | calder_random_groups | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.8 | calder_random_groups | calder | auc_hits | 24.3 | -7.902 | 61.9 | 10 |
| synthetic | A | 0.8 | calder_uninformative_features | calder | confirmed_hits | -14.9 | -17.6 | -12.4 | 10 |
| synthetic | A | 0.8 | calder_uninformative_features | calder | group_coverage | 0.1 | -0.4 | 0.6 | 10 |
| synthetic | A | 0.8 | calder_uninformative_features | calder | auc_hits | -307.8 | -417.5 | -203.3 | 10 |
| synthetic | A | 0.9 | calder_no_calibration | calder | confirmed_hits | -1.8 | -3.1 | -0.2 | 10 |
| synthetic | A | 0.9 | calder_no_calibration | calder | group_coverage | 0 | -0.3 | 0.3 | 10 |
| synthetic | A | 0.9 | calder_no_calibration | calder | auc_hits | -64.7 | -108.3 | -13.89 | 10 |
| synthetic | A | 0.9 | calder_no_diversity | calder | confirmed_hits | 0.9 | -0.8 | 2.6 | 10 |
| synthetic | A | 0.9 | calder_no_diversity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.9 | calder_no_diversity | calder | auc_hits | 46 | 15.7 | 80.3 | 10 |
| synthetic | A | 0.9 | calder_no_similarity | calder | confirmed_hits | 1.5 | -0.4 | 3.1 | 10 |
| synthetic | A | 0.9 | calder_no_similarity | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.9 | calder_no_similarity | calder | auc_hits | 72.7 | 39.6 | 100.5 | 10 |
| synthetic | A | 0.9 | calder_no_uncertainty | calder | confirmed_hits | -0.1 | -1.7 | 1.502 | 10 |
| synthetic | A | 0.9 | calder_no_uncertainty | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | A | 0.9 | calder_no_uncertainty | calder | auc_hits | 30.2 | 8.5 | 51.2 | 10 |
| synthetic | A | 0.9 | calder_random_groups | calder | confirmed_hits | 1.5 | -0.2 | 3.1 | 10 |
| synthetic | A | 0.9 | calder_random_groups | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | A | 0.9 | calder_random_groups | calder | auc_hits | 63.5 | 36.9 | 93.4 | 10 |
| synthetic | A | 0.9 | calder_uninformative_features | calder | confirmed_hits | -9.2 | -11.4 | -7.1 | 10 |
| synthetic | A | 0.9 | calder_uninformative_features | calder | group_coverage | -0.2 | -0.6 | 0.2 | 10 |
| synthetic | A | 0.9 | calder_uninformative_features | calder | auc_hits | -197.5 | -276.3 | -124.2 | 10 |
| synthetic | B | 0.7 | calder_no_calibration | calder | confirmed_hits | 1.4 | -0.4 | 3.1 | 10 |
| synthetic | B | 0.7 | calder_no_calibration | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.7 | calder_no_calibration | calder | auc_hits | 19.1 | -19.5 | 60.6 | 10 |
| synthetic | B | 0.7 | calder_no_diversity | calder | confirmed_hits | 0 | 0 | 0 | 10 |
| synthetic | B | 0.7 | calder_no_diversity | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | B | 0.7 | calder_no_diversity | calder | auc_hits | 0.6 | -5.1 | 6.9 | 10 |
| synthetic | B | 0.7 | calder_no_similarity | calder | confirmed_hits | 0.9 | -0.8 | 2.5 | 10 |
| synthetic | B | 0.7 | calder_no_similarity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.7 | calder_no_similarity | calder | auc_hits | 54.3 | 15.3 | 91.6 | 10 |
| synthetic | B | 0.7 | calder_no_uncertainty | calder | confirmed_hits | -0.2 | -1.7 | 1.3 | 10 |
| synthetic | B | 0.7 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.7 | calder_no_uncertainty | calder | auc_hits | -2.9 | -31.21 | 26.4 | 10 |
| synthetic | B | 0.7 | calder_random_groups | calder | confirmed_hits | -0.1 | -0.9 | 0.6 | 10 |
| synthetic | B | 0.7 | calder_random_groups | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | B | 0.7 | calder_random_groups | calder | auc_hits | 5.8 | -8.6 | 28.6 | 10 |
| synthetic | B | 0.7 | calder_uninformative_features | calder | confirmed_hits | 1.3 | -0.6 | 3.4 | 10 |
| synthetic | B | 0.7 | calder_uninformative_features | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.7 | calder_uninformative_features | calder | auc_hits | 43.7 | -24.13 | 100.9 | 10 |
| synthetic | B | 0.8 | calder_no_calibration | calder | confirmed_hits | 1.1 | -1.5 | 3.7 | 10 |
| synthetic | B | 0.8 | calder_no_calibration | calder | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | B | 0.8 | calder_no_calibration | calder | auc_hits | 17.6 | -16.6 | 62.8 | 10 |
| synthetic | B | 0.8 | calder_no_diversity | calder | confirmed_hits | -0.3 | -0.9 | 0.2 | 10 |
| synthetic | B | 0.8 | calder_no_diversity | calder | group_coverage | -0.4 | -0.8025 | 0 | 10 |
| synthetic | B | 0.8 | calder_no_diversity | calder | auc_hits | -7.1 | -17 | 0 | 10 |
| synthetic | B | 0.8 | calder_no_similarity | calder | confirmed_hits | 1.3 | -0.4 | 3 | 10 |
| synthetic | B | 0.8 | calder_no_similarity | calder | group_coverage | -0.3 | -0.7025 | 0 | 10 |
| synthetic | B | 0.8 | calder_no_similarity | calder | auc_hits | 20.6 | -16.2 | 56.2 | 10 |
| synthetic | B | 0.8 | calder_no_uncertainty | calder | confirmed_hits | 0.2 | -0.8 | 1.1 | 10 |
| synthetic | B | 0.8 | calder_no_uncertainty | calder | group_coverage | -0.3 | -0.8 | 0 | 10 |
| synthetic | B | 0.8 | calder_no_uncertainty | calder | auc_hits | 5.5 | -8.2 | 21.1 | 10 |
| synthetic | B | 0.8 | calder_random_groups | calder | confirmed_hits | 0.2 | -0.5 | 1.1 | 10 |
| synthetic | B | 0.8 | calder_random_groups | calder | group_coverage | -0.3 | -0.7025 | 0 | 10 |
| synthetic | B | 0.8 | calder_random_groups | calder | auc_hits | -7.6 | -22.71 | 6.6 | 10 |
| synthetic | B | 0.8 | calder_uninformative_features | calder | confirmed_hits | -1.4 | -2.8 | 0 | 10 |
| synthetic | B | 0.8 | calder_uninformative_features | calder | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | B | 0.8 | calder_uninformative_features | calder | auc_hits | -39.9 | -81.8 | -0.3975 | 10 |
| synthetic | B | 0.9 | calder_no_calibration | calder | confirmed_hits | 0.2 | -0.3 | 0.8 | 10 |
| synthetic | B | 0.9 | calder_no_calibration | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | B | 0.9 | calder_no_calibration | calder | auc_hits | 0.2 | -12.9 | 13.1 | 10 |
| synthetic | B | 0.9 | calder_no_diversity | calder | confirmed_hits | 0.5 | -0.2 | 1.1 | 10 |
| synthetic | B | 0.9 | calder_no_diversity | calder | group_coverage | -0.2 | -0.6 | 0.2 | 10 |
| synthetic | B | 0.9 | calder_no_diversity | calder | auc_hits | 22.6 | 1.498 | 50.5 | 10 |
| synthetic | B | 0.9 | calder_no_similarity | calder | confirmed_hits | 1.3 | 0 | 2.6 | 10 |
| synthetic | B | 0.9 | calder_no_similarity | calder | group_coverage | 0.3 | 0.1 | 0.6 | 10 |
| synthetic | B | 0.9 | calder_no_similarity | calder | auc_hits | 11 | -18.4 | 43.2 | 10 |
| synthetic | B | 0.9 | calder_no_uncertainty | calder | confirmed_hits | -0.2 | -1.2 | 0.6 | 10 |
| synthetic | B | 0.9 | calder_no_uncertainty | calder | group_coverage | 0 | -0.4 | 0.4 | 10 |
| synthetic | B | 0.9 | calder_no_uncertainty | calder | auc_hits | 2.4 | -27.6 | 31.7 | 10 |
| synthetic | B | 0.9 | calder_random_groups | calder | confirmed_hits | -0.3 | -1.602 | 1 | 10 |
| synthetic | B | 0.9 | calder_random_groups | calder | group_coverage | -0.1 | -0.5 | 0.3 | 10 |
| synthetic | B | 0.9 | calder_random_groups | calder | auc_hits | -3.5 | -32.5 | 24.7 | 10 |
| synthetic | B | 0.9 | calder_uninformative_features | calder | confirmed_hits | -1.4 | -3.402 | 0.4 | 10 |
| synthetic | B | 0.9 | calder_uninformative_features | calder | group_coverage | 0.1 | -0.6 | 0.8 | 10 |
| synthetic | B | 0.9 | calder_uninformative_features | calder | auc_hits | -42.5 | -100.7 | 9.007 | 10 |
| synthetic | C | 0.7 | calder_no_calibration | calder | confirmed_hits | 1.4 | -0.1 | 2.9 | 10 |
| synthetic | C | 0.7 | calder_no_calibration | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_no_calibration | calder | auc_hits | 35.3 | 11.7 | 59 | 10 |
| synthetic | C | 0.7 | calder_no_diversity | calder | confirmed_hits | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.7 | calder_no_diversity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_no_diversity | calder | auc_hits | -2.3 | -6.9 | 0 | 10 |
| synthetic | C | 0.7 | calder_no_similarity | calder | confirmed_hits | 0.1 | -1.4 | 1.6 | 10 |
| synthetic | C | 0.7 | calder_no_similarity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_no_similarity | calder | auc_hits | 7.9 | -37.1 | 52.31 | 10 |
| synthetic | C | 0.7 | calder_no_uncertainty | calder | confirmed_hits | 0.7 | -1.1 | 2.5 | 10 |
| synthetic | C | 0.7 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_no_uncertainty | calder | auc_hits | 25.8 | -25.6 | 80.6 | 10 |
| synthetic | C | 0.7 | calder_random_groups | calder | confirmed_hits | 0.2 | 0 | 0.5 | 10 |
| synthetic | C | 0.7 | calder_random_groups | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_random_groups | calder | auc_hits | 5.4 | -6.9 | 23.1 | 10 |
| synthetic | C | 0.7 | calder_uninformative_features | calder | confirmed_hits | -1.3 | -3.4 | 1 | 10 |
| synthetic | C | 0.7 | calder_uninformative_features | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.7 | calder_uninformative_features | calder | auc_hits | -41.5 | -109.6 | 24.51 | 10 |
| synthetic | C | 0.8 | calder_no_calibration | calder | confirmed_hits | 2.1 | 0.6 | 3.6 | 10 |
| synthetic | C | 0.8 | calder_no_calibration | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.8 | calder_no_calibration | calder | auc_hits | 31.1 | 11.1 | 51.3 | 10 |
| synthetic | C | 0.8 | calder_no_diversity | calder | confirmed_hits | 0.4 | 0 | 1.2 | 10 |
| synthetic | C | 0.8 | calder_no_diversity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | calder_no_diversity | calder | auc_hits | 2.6 | -3.3 | 11.1 | 10 |
| synthetic | C | 0.8 | calder_no_similarity | calder | confirmed_hits | 0.9 | -1.1 | 2.7 | 10 |
| synthetic | C | 0.8 | calder_no_similarity | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.8 | calder_no_similarity | calder | auc_hits | 36.4 | 7.598 | 67.2 | 10 |
| synthetic | C | 0.8 | calder_no_uncertainty | calder | confirmed_hits | 1 | -0.9 | 2.8 | 10 |
| synthetic | C | 0.8 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | calder_no_uncertainty | calder | auc_hits | 42 | 8.798 | 72.4 | 10 |
| synthetic | C | 0.8 | calder_random_groups | calder | confirmed_hits | 0.1 | -0.9 | 1.2 | 10 |
| synthetic | C | 0.8 | calder_random_groups | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.8 | calder_random_groups | calder | auc_hits | 1.8 | -4.1 | 10.3 | 10 |
| synthetic | C | 0.8 | calder_uninformative_features | calder | confirmed_hits | 0.4 | -1.5 | 2.4 | 10 |
| synthetic | C | 0.8 | calder_uninformative_features | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.8 | calder_uninformative_features | calder | auc_hits | -4.9 | -52.01 | 36.8 | 10 |
| synthetic | C | 0.9 | calder_no_calibration | calder | confirmed_hits | 0.9 | -0.5 | 2.1 | 10 |
| synthetic | C | 0.9 | calder_no_calibration | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.9 | calder_no_calibration | calder | auc_hits | 4.2 | -31.3 | 39.1 | 10 |
| synthetic | C | 0.9 | calder_no_diversity | calder | confirmed_hits | 0.4 | -0.4 | 1.005 | 10 |
| synthetic | C | 0.9 | calder_no_diversity | calder | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | C | 0.9 | calder_no_diversity | calder | auc_hits | 1.5 | -10.3 | 14.7 | 10 |
| synthetic | C | 0.9 | calder_no_similarity | calder | confirmed_hits | 1.1 | -0.9 | 2.8 | 10 |
| synthetic | C | 0.9 | calder_no_similarity | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.9 | calder_no_similarity | calder | auc_hits | 26 | -3.705 | 60.5 | 10 |
| synthetic | C | 0.9 | calder_no_uncertainty | calder | confirmed_hits | 0.5 | -0.5025 | 1.7 | 10 |
| synthetic | C | 0.9 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | C | 0.9 | calder_no_uncertainty | calder | auc_hits | 14.2 | -13 | 45.2 | 10 |
| synthetic | C | 0.9 | calder_random_groups | calder | confirmed_hits | 0.5 | -0.2 | 1.2 | 10 |
| synthetic | C | 0.9 | calder_random_groups | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | C | 0.9 | calder_random_groups | calder | auc_hits | 15.1 | 1.7 | 28.6 | 10 |
| synthetic | C | 0.9 | calder_uninformative_features | calder | confirmed_hits | -0.4 | -2.5 | 1.3 | 10 |
| synthetic | C | 0.9 | calder_uninformative_features | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | C | 0.9 | calder_uninformative_features | calder | auc_hits | -31.9 | -88 | 20.1 | 10 |
| synthetic | D | 0.7 | calder_no_calibration | calder | confirmed_hits | -0.3 | -3.1 | 2.2 | 10 |
| synthetic | D | 0.7 | calder_no_calibration | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_no_calibration | calder | auc_hits | -1.4 | -46.1 | 41.8 | 10 |
| synthetic | D | 0.7 | calder_no_diversity | calder | confirmed_hits | -0.2 | -1.1 | 0.7 | 10 |
| synthetic | D | 0.7 | calder_no_diversity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_no_diversity | calder | auc_hits | -5.8 | -16.7 | 1.4 | 10 |
| synthetic | D | 0.7 | calder_no_similarity | calder | confirmed_hits | 0.5 | -1.7 | 2.8 | 10 |
| synthetic | D | 0.7 | calder_no_similarity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_no_similarity | calder | auc_hits | -5.6 | -75.81 | 56.52 | 10 |
| synthetic | D | 0.7 | calder_no_uncertainty | calder | confirmed_hits | -0.1 | -1.9 | 2 | 10 |
| synthetic | D | 0.7 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_no_uncertainty | calder | auc_hits | -5.7 | -45.9 | 50.42 | 10 |
| synthetic | D | 0.7 | calder_random_groups | calder | confirmed_hits | -0.2 | -1.1 | 0.7 | 10 |
| synthetic | D | 0.7 | calder_random_groups | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_random_groups | calder | auc_hits | -5.8 | -16.7 | 1.4 | 10 |
| synthetic | D | 0.7 | calder_uninformative_features | calder | confirmed_hits | -0.4 | -3.8 | 2.8 | 10 |
| synthetic | D | 0.7 | calder_uninformative_features | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.7 | calder_uninformative_features | calder | auc_hits | -28.8 | -104.5 | 36.41 | 10 |
| synthetic | D | 0.8 | calder_no_calibration | calder | confirmed_hits | 0.7 | -1.2 | 2.1 | 10 |
| synthetic | D | 0.8 | calder_no_calibration | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | calder_no_calibration | calder | auc_hits | 16.4 | -17.6 | 56.5 | 10 |
| synthetic | D | 0.8 | calder_no_diversity | calder | confirmed_hits | -0.8 | -2.2 | 0.4 | 10 |
| synthetic | D | 0.8 | calder_no_diversity | calder | group_coverage | -0.1 | -0.3 | 0 | 10 |
| synthetic | D | 0.8 | calder_no_diversity | calder | auc_hits | -21.2 | -44.7 | -2.998 | 10 |
| synthetic | D | 0.8 | calder_no_similarity | calder | confirmed_hits | 0.3 | -1.6 | 2.2 | 10 |
| synthetic | D | 0.8 | calder_no_similarity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | calder_no_similarity | calder | auc_hits | 2.9 | -34.3 | 41.7 | 10 |
| synthetic | D | 0.8 | calder_no_uncertainty | calder | confirmed_hits | 0.4 | -0.6 | 1.4 | 10 |
| synthetic | D | 0.8 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | calder_no_uncertainty | calder | auc_hits | 11.8 | -10.2 | 34.8 | 10 |
| synthetic | D | 0.8 | calder_random_groups | calder | confirmed_hits | -0.7 | -2.7 | 1 | 10 |
| synthetic | D | 0.8 | calder_random_groups | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.8 | calder_random_groups | calder | auc_hits | -14.5 | -41.5 | 8.4 | 10 |
| synthetic | D | 0.8 | calder_uninformative_features | calder | confirmed_hits | -1.1 | -3.8 | 1.502 | 10 |
| synthetic | D | 0.8 | calder_uninformative_features | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.8 | calder_uninformative_features | calder | auc_hits | -41.9 | -116.6 | 22.9 | 10 |
| synthetic | D | 0.9 | calder_no_calibration | calder | confirmed_hits | 0.8 | -0.7 | 2.4 | 10 |
| synthetic | D | 0.9 | calder_no_calibration | calder | group_coverage | 0.1 | 0 | 0.3 | 10 |
| synthetic | D | 0.9 | calder_no_calibration | calder | auc_hits | 9.5 | -33.01 | 56.6 | 10 |
| synthetic | D | 0.9 | calder_no_diversity | calder | confirmed_hits | -0.3 | -1 | 0.3025 | 10 |
| synthetic | D | 0.9 | calder_no_diversity | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.9 | calder_no_diversity | calder | auc_hits | -26.1 | -47 | -8.4 | 10 |
| synthetic | D | 0.9 | calder_no_similarity | calder | confirmed_hits | 0.8 | -0.6 | 2.2 | 10 |
| synthetic | D | 0.9 | calder_no_similarity | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.9 | calder_no_similarity | calder | auc_hits | 0.3 | -25.9 | 25.41 | 10 |
| synthetic | D | 0.9 | calder_no_uncertainty | calder | confirmed_hits | 0 | -0.8 | 0.8 | 10 |
| synthetic | D | 0.9 | calder_no_uncertainty | calder | group_coverage | 0 | 0 | 0 | 10 |
| synthetic | D | 0.9 | calder_no_uncertainty | calder | auc_hits | 3.8 | -23.2 | 34.6 | 10 |
| synthetic | D | 0.9 | calder_random_groups | calder | confirmed_hits | -0.9 | -1.8 | -0.1 | 10 |
| synthetic | D | 0.9 | calder_random_groups | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.9 | calder_random_groups | calder | auc_hits | -22 | -40.6 | -3.6 | 10 |
| synthetic | D | 0.9 | calder_uninformative_features | calder | confirmed_hits | -0.6 | -2.2 | 1.1 | 10 |
| synthetic | D | 0.9 | calder_uninformative_features | calder | group_coverage | -0.2 | -0.5 | 0 | 10 |
| synthetic | D | 0.9 | calder_uninformative_features | calder | auc_hits | -26.2 | -58.6 | 13 | 10 |

Full confidence intervals are in `results.csv`; seed endpoints are in `seed_level_results.csv`. Prequential latent/measured diagnostics remain separate in `calibration.csv`. Figure availability and links are recorded in [summary.md](summary.md).

## 7. Ablation studies

Completed controlled variants are preserved in `ablations.csv`; component differences must be assessed within each fixed scenario and threshold.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean |
| --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.7 | calder | 10 | 43 | 4.4 |
| synthetic | A | 0.7 | calder_no_calibration | 10 | 43.4 | 4.2 |
| synthetic | A | 0.7 | calder_no_diversity | 10 | 44.4 | 3.7 |
| synthetic | A | 0.7 | calder_no_similarity | 10 | 45.9 | 4.1 |
| synthetic | A | 0.7 | calder_no_uncertainty | 10 | 41.6 | 4.1 |
| synthetic | A | 0.7 | calder_random_groups | 10 | 44.4 | 3.7 |
| synthetic | A | 0.7 | calder_uninformative_features | 10 | 29.5 | 4.6 |
| synthetic | A | 0.8 | calder | 10 | 32.6 | 3.1 |
| synthetic | A | 0.8 | calder_no_calibration | 10 | 32.1 | 3.2 |
| synthetic | A | 0.8 | calder_no_diversity | 10 | 33.2 | 3.2 |
| synthetic | A | 0.8 | calder_no_similarity | 10 | 36.3 | 3.2 |
| synthetic | A | 0.8 | calder_no_uncertainty | 10 | 32.5 | 3.1 |
| synthetic | A | 0.8 | calder_random_groups | 10 | 32.8 | 3.2 |
| synthetic | A | 0.8 | calder_uninformative_features | 10 | 17.7 | 3.2 |
| synthetic | A | 0.9 | calder | 10 | 18.7 | 1.6 |
| synthetic | A | 0.9 | calder_no_calibration | 10 | 16.9 | 1.6 |
| synthetic | A | 0.9 | calder_no_diversity | 10 | 19.6 | 1.6 |
| synthetic | A | 0.9 | calder_no_similarity | 10 | 20.2 | 1.7 |
| synthetic | A | 0.9 | calder_no_uncertainty | 10 | 18.6 | 1.7 |
| synthetic | A | 0.9 | calder_random_groups | 10 | 20.2 | 1.6 |
| synthetic | A | 0.9 | calder_uninformative_features | 10 | 9.5 | 1.4 |
| synthetic | B | 0.7 | calder | 10 | 31.5 | 5 |
| synthetic | B | 0.7 | calder_no_calibration | 10 | 32.9 | 5 |
| synthetic | B | 0.7 | calder_no_diversity | 10 | 31.5 | 4.8 |
| synthetic | B | 0.7 | calder_no_similarity | 10 | 32.4 | 5 |
| synthetic | B | 0.7 | calder_no_uncertainty | 10 | 31.3 | 5 |
| synthetic | B | 0.7 | calder_random_groups | 10 | 31.4 | 4.8 |
| synthetic | B | 0.7 | calder_uninformative_features | 10 | 32.8 | 5 |
| synthetic | B | 0.8 | calder | 10 | 20.5 | 5 |
| synthetic | B | 0.8 | calder_no_calibration | 10 | 21.6 | 4.9 |
| synthetic | B | 0.8 | calder_no_diversity | 10 | 20.2 | 4.6 |
| synthetic | B | 0.8 | calder_no_similarity | 10 | 21.8 | 4.7 |
| synthetic | B | 0.8 | calder_no_uncertainty | 10 | 20.7 | 4.7 |
| synthetic | B | 0.8 | calder_random_groups | 10 | 20.7 | 4.7 |
| synthetic | B | 0.8 | calder_uninformative_features | 10 | 19.1 | 4.9 |
| synthetic | B | 0.9 | calder | 10 | 10.9 | 3.9 |
| synthetic | B | 0.9 | calder_no_calibration | 10 | 11.1 | 3.9 |
| synthetic | B | 0.9 | calder_no_diversity | 10 | 11.4 | 3.7 |
| synthetic | B | 0.9 | calder_no_similarity | 10 | 12.2 | 4.2 |
| synthetic | B | 0.9 | calder_no_uncertainty | 10 | 10.7 | 3.9 |
| synthetic | B | 0.9 | calder_random_groups | 10 | 10.6 | 3.8 |
| synthetic | B | 0.9 | calder_uninformative_features | 10 | 9.5 | 4 |
| synthetic | C | 0.7 | calder | 10 | 33 | 5 |
| synthetic | C | 0.7 | calder_no_calibration | 10 | 34.4 | 5 |
| synthetic | C | 0.7 | calder_no_diversity | 10 | 33.1 | 5 |
| synthetic | C | 0.7 | calder_no_similarity | 10 | 33.1 | 5 |
| synthetic | C | 0.7 | calder_no_uncertainty | 10 | 33.7 | 5 |
| synthetic | C | 0.7 | calder_random_groups | 10 | 33.2 | 5 |
| synthetic | C | 0.7 | calder_uninformative_features | 10 | 31.7 | 5 |
| synthetic | C | 0.8 | calder | 10 | 23.6 | 4.9 |
| synthetic | C | 0.8 | calder_no_calibration | 10 | 25.7 | 5 |
| synthetic | C | 0.8 | calder_no_diversity | 10 | 24 | 4.9 |
| synthetic | C | 0.8 | calder_no_similarity | 10 | 24.5 | 5 |
| synthetic | C | 0.8 | calder_no_uncertainty | 10 | 24.6 | 4.9 |
| synthetic | C | 0.8 | calder_random_groups | 10 | 23.7 | 4.9 |
| synthetic | C | 0.8 | calder_uninformative_features | 10 | 24 | 5 |
| synthetic | C | 0.9 | calder | 10 | 16.1 | 4.9 |
| synthetic | C | 0.9 | calder_no_calibration | 10 | 17 | 5 |
| synthetic | C | 0.9 | calder_no_diversity | 10 | 16.5 | 4.8 |
| synthetic | C | 0.9 | calder_no_similarity | 10 | 17.2 | 5 |
| synthetic | C | 0.9 | calder_no_uncertainty | 10 | 16.6 | 4.9 |
| synthetic | C | 0.9 | calder_random_groups | 10 | 16.6 | 4.7 |
| synthetic | C | 0.9 | calder_uninformative_features | 10 | 15.7 | 5 |
| synthetic | D | 0.7 | calder | 10 | 32.1 | 5 |
| synthetic | D | 0.7 | calder_no_calibration | 10 | 31.8 | 5 |
| synthetic | D | 0.7 | calder_no_diversity | 10 | 31.9 | 5 |
| synthetic | D | 0.7 | calder_no_similarity | 10 | 32.6 | 5 |
| synthetic | D | 0.7 | calder_no_uncertainty | 10 | 32 | 5 |
| synthetic | D | 0.7 | calder_random_groups | 10 | 31.9 | 5 |
| synthetic | D | 0.7 | calder_uninformative_features | 10 | 31.7 | 5 |
| synthetic | D | 0.8 | calder | 10 | 22.5 | 5 |
| synthetic | D | 0.8 | calder_no_calibration | 10 | 23.2 | 5 |
| synthetic | D | 0.8 | calder_no_diversity | 10 | 21.7 | 4.9 |
| synthetic | D | 0.8 | calder_no_similarity | 10 | 22.8 | 5 |
| synthetic | D | 0.8 | calder_no_uncertainty | 10 | 22.9 | 5 |
| synthetic | D | 0.8 | calder_random_groups | 10 | 21.8 | 4.8 |
| synthetic | D | 0.8 | calder_uninformative_features | 10 | 21.4 | 5 |
| synthetic | D | 0.9 | calder | 10 | 13.1 | 4.9 |
| synthetic | D | 0.9 | calder_no_calibration | 10 | 13.9 | 5 |
| synthetic | D | 0.9 | calder_no_diversity | 10 | 12.8 | 4.7 |
| synthetic | D | 0.9 | calder_no_similarity | 10 | 13.9 | 4.9 |
| synthetic | D | 0.9 | calder_no_uncertainty | 10 | 13.1 | 4.9 |
| synthetic | D | 0.9 | calder_random_groups | 10 | 12.2 | 4.7 |
| synthetic | D | 0.9 | calder_uninformative_features | 10 | 12.5 | 4.7 |

## 8. Limitations

Noisy pseudobulk response magnitudes depend on filtering, feature choice and cell sampling. Bootstrap uncertainty omits some experimental confounding. Shared controls can correlate estimation errors. Retrospective replay provides no prospective laboratory validation. Prequential examples are actively selected and therefore biased toward the acquisition distribution. Small pilot calibration reserves limit calibration power and can force fallback. Interval coverage is descriptive and has no exchangeability-based guarantee under adaptive selection.

Network annotations can be incomplete, inaccurate or context-insensitive. Feature clusters do not establish biological mechanisms. Conditional independence in coverage utility is approximate. RPE1 reruns evaluate a new context, whereas direct K562-to-RPE1 predictive transfer is a different optional experiment. Synthetic success cannot establish biological performance. A high response score does not establish therapeutic usefulness, and intervention effects do not automatically establish direct causal regulatory edges.

## 9. Discussion

Interpret hit discovery and group coverage jointly without combining them into an arbitrary headline score. If diversity adds groups while reducing hits, quantify both changes. If calibration improves Brier score without improving acquisition, examine monotone ranking preservation, sample size, shifted score distributions and saturation of the coverage term. A more accurate model may still select fewer hits under a particular finite-budget rule.

## 10. Conclusion

We evaluated budgeted perturbation selection using 840 saved trajectories, 7 selection/reference rules and 10 random seeds in controlled synthetic simulation. Seed-paired bootstrap intervals quantify algorithmic variation conditional on the fixed pools. No laboratory experiments were performed, and these results do not establish universal policy superiority. Biological K562–RPE1 performance remains unmeasured in this run.

## 11. References

See the numbered [bibliography](bibliography.md).
