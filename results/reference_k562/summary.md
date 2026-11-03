# CALDER experiment summary

This report was generated from saved experiment outputs. No real laboratory experiments were performed by CALDER.

Experiment type: retrospective replay of public perturbation measurements.

Configuration hash: `f22d66d39fdb00d01c7278aefdc90bc7bd0c5d229826df854e0b7227b475e1b1`. Git commit: `7e95d257bd7aaec8247dae6ff8f8d524fdb2c024`.

## Measured endpoint results

Confirmed hits use observed log effects at the pilot-derived threshold. Synthetic latent-hit counts, when present, are a distinct ground-truth diagnostic. Group coverage counts predefined groups containing at least one confirmed measured hit.

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

Full bootstrap intervals and secondary metrics are in `results.csv`; every run endpoint is in `seed_level_results.csv`. Seeds describe algorithmic variation on a fixed experimental pool, not independent biological replications.

## Paired policy differences

Positive differences favour the named policy. Intervals resample whole paired seeds within each dataset, scenario and threshold. They quantify Monte Carlo variation conditional on this pool. Multiple contrasts are descriptive; no multiplicity-adjusted superiority claim is made.

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

## Calibration and prediction

`calibration.csv` reports raw, pilot-calibrated and acquisition-used probability diagnostics: Brier score, fixed-bin ECE, prediction error, uncertainty–error association and interval coverage. Bin counts and group-specific diagnostics are in `reliability.csv` and `grouped_calibration.csv`. These are pre-query predictions on adaptively acquired examples. Selection bias prevents interpreting them as exchangeable population-level validation. Latent probabilities evaluated against noisy measured hit labels also reflect label noise. Calibration fitted to measured pilot hit labels is a measured-label correction, not guaranteed latent-effect calibration.

Measured predictive intervals add the recorded sampling variance and configured likelihood-discrepancy variance floor to epistemic variance. Their target is labelled `measured_predictive_model_assumed`; this is a model assumption and is heuristic for ridge or ensemble uncertainty. Older saved runs restore the floor from their resolved configuration in memory, preserving original prediction tables. Latent-effect intervals use epistemic variance alone.

## Figure availability

- figure_1: generated. [view](figures/figure_1_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.7.png), [view](figures/figure_1_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png), [view](figures/figure_1_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.9.png)

- figure_2: generated. [view](figures/figure_2_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.7.png), [view](figures/figure_2_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png), [view](figures/figure_2_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.9.png)

- figure_3: generated. [view](figures/figure_3_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.7.png), [view](figures/figure_3_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png), [view](figures/figure_3_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.9.png)

- figure_4: generated. [view](figures/figure_4_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.7.png), [view](figures/figure_4_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png), [view](figures/figure_4_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.9.png)

- figure_5: unavailable. Required saved inputs absent

- figure_6: unavailable. Both K562 and RPE1 measured benchmark endpoints are required

- figure_7: generated. [view](figures/figure_7_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.7.png), [view](figures/figure_7_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.8.png), [view](figures/figure_7_dataset-replogle22k562__scenario-replogle22k562__threshold_quantile-0.9.png)

- figure_8: generated. [view](figures/figure_8_dataset-replogle22k562__scenario-replogle22k562.png)


## Interpretation limits

The policy is not guaranteed to outperform Probability-of-Hit or other baselines. Negative and mixed outcomes remain in the saved tables. A high transcriptional-response score does not establish therapeutic usefulness, and perturbation-induced expression changes do not identify direct causal regulatory edges. STRING communities are prior association groups; control-feature clusters are artificial feature clusters. Projected feature geometry has no direct causal interpretation.

K562–RPE1 comparison requires measured outputs from both contexts under frozen method hyperparameters; missing context results are not replaced with synthetic values. Offline oracle rankings are full-information upper-bound references and are not deployable policies.


Sources and protocol are documented in the repository's `reports/bibliography.md` and `reports/research_report.md`.
