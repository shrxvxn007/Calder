# CALDER experiment summary

This report was generated from saved experiment outputs. No real laboratory experiments were performed by CALDER.

Experiment type: controlled synthetic simulation.

Configuration hash: `0da38cfefe2b730481b801499fb1140bbbff77ffe18958d70cee4ec2fa36d436`. Git commit: `15de4bebfb212eca645a74806bafbee35beb3a5b`.

## Measured endpoint results

Confirmed hits use observed log effects at the pilot-derived threshold. Synthetic latent-hit counts, when present, are a distinct ground-truth diagnostic. Group coverage counts predefined groups containing at least one confirmed measured hit.

| dataset | scenario | threshold_quantile | policy | n_seeds | confirmed_hits_mean | group_coverage_mean | auc_hits_mean |
| --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.8 | calder | 3 | 32.33 | 2.667 | 958.3 |
| synthetic | A | 0.8 | probability_of_hit | 3 | 37.67 | 2.667 | 1228 |
| synthetic | A | 0.8 | random | 3 | 22.67 | 3 | 810.3 |
| synthetic | B | 0.8 | calder | 3 | 21.33 | 5 | 748 |
| synthetic | B | 0.8 | probability_of_hit | 3 | 24.67 | 4.667 | 792.3 |
| synthetic | B | 0.8 | random | 3 | 20 | 5 | 740.3 |
| synthetic | C | 0.8 | calder | 3 | 24 | 5 | 806.7 |
| synthetic | C | 0.8 | probability_of_hit | 3 | 25.67 | 5 | 856.3 |
| synthetic | C | 0.8 | random | 3 | 24.67 | 5 | 860.3 |
| synthetic | D | 0.8 | calder | 3 | 23 | 5 | 789.7 |
| synthetic | D | 0.8 | probability_of_hit | 3 | 23.67 | 5 | 835 |
| synthetic | D | 0.8 | random | 3 | 19 | 5 | 694.3 |

Full bootstrap intervals and secondary metrics are in `results.csv`; every run endpoint is in `seed_level_results.csv`. Seeds describe algorithmic variation on a fixed experimental pool, not independent biological replications.

## Paired policy differences

Positive differences favour the named policy. Intervals resample whole paired seeds within each dataset, scenario and threshold. They quantify Monte Carlo variation conditional on this pool. Multiple contrasts are descriptive; no multiplicity-adjusted superiority claim is made.

| dataset | scenario | threshold_quantile | policy | reference | metric | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| synthetic | A | 0.8 | calder | probability_of_hit | confirmed_hits | -5.333 | -9 | -3 | 3 |
| synthetic | A | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | A | 0.8 | calder | probability_of_hit | auc_hits | -269.7 | -342 | -172 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | confirmed_hits | -15 | -21 | -10 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | A | 0.8 | random | probability_of_hit | auc_hits | -417.7 | -572 | -310 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | confirmed_hits | -3.333 | -4 | -2 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | B | 0.8 | calder | probability_of_hit | auc_hits | -44.33 | -118 | 2 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | confirmed_hits | -4.667 | -7 | -3 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | group_coverage | 0.3333 | 0 | 1 | 3 |
| synthetic | B | 0.8 | random | probability_of_hit | auc_hits | -52 | -108 | 20 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | confirmed_hits | -1.667 | -3 | -1 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | C | 0.8 | calder | probability_of_hit | auc_hits | -49.67 | -109 | 18 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | confirmed_hits | -1 | -2 | 0 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | C | 0.8 | random | probability_of_hit | auc_hits | 4 | -106 | 95 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | confirmed_hits | -0.6667 | -2 | 1 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | D | 0.8 | calder | probability_of_hit | auc_hits | -45.33 | -129 | 30 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | confirmed_hits | -4.667 | -8 | -3 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | group_coverage | 0 | 0 | 0 | 3 |
| synthetic | D | 0.8 | random | probability_of_hit | auc_hits | -140.7 | -211 | -25 | 3 |

## Calibration and prediction

`calibration.csv` reports raw, pilot-calibrated and acquisition-used probability diagnostics: Brier score, fixed-bin ECE, prediction error, uncertainty–error association and interval coverage. Bin counts and group-specific diagnostics are in `reliability.csv` and `grouped_calibration.csv`. These are pre-query predictions on adaptively acquired examples. Selection bias prevents interpreting them as exchangeable population-level validation. Latent probabilities evaluated against noisy measured hit labels also reflect label noise. Calibration fitted to measured pilot hit labels is a measured-label correction, not guaranteed latent-effect calibration.

## Figure availability

- figure_1: generated. [view](figures/figure_1_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [view](figures/figure_1_dataset-synthetic__scenario-B__threshold_quantile-0.8.png), [view](figures/figure_1_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [view](figures/figure_1_dataset-synthetic__scenario-D__threshold_quantile-0.8.png)

- figure_2: generated. [view](figures/figure_2_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [view](figures/figure_2_dataset-synthetic__scenario-B__threshold_quantile-0.8.png), [view](figures/figure_2_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [view](figures/figure_2_dataset-synthetic__scenario-D__threshold_quantile-0.8.png)

- figure_3: generated. [view](figures/figure_3_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [view](figures/figure_3_dataset-synthetic__scenario-B__threshold_quantile-0.8.png), [view](figures/figure_3_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [view](figures/figure_3_dataset-synthetic__scenario-D__threshold_quantile-0.8.png)

- figure_4: generated. [view](figures/figure_4_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [view](figures/figure_4_dataset-synthetic__scenario-B__threshold_quantile-0.8.png), [view](figures/figure_4_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [view](figures/figure_4_dataset-synthetic__scenario-D__threshold_quantile-0.8.png)

- figure_5: unavailable. Required saved inputs absent

- figure_6: unavailable. Both K562 and RPE1 measured benchmark endpoints are required

- figure_7: generated. [view](figures/figure_7_dataset-synthetic__scenario-A__threshold_quantile-0.8.png), [view](figures/figure_7_dataset-synthetic__scenario-B__threshold_quantile-0.8.png), [view](figures/figure_7_dataset-synthetic__scenario-C__threshold_quantile-0.8.png), [view](figures/figure_7_dataset-synthetic__scenario-D__threshold_quantile-0.8.png)

- figure_8: generated. [view](figures/figure_8_dataset-synthetic__scenario-A.png), [view](figures/figure_8_dataset-synthetic__scenario-B.png), [view](figures/figure_8_dataset-synthetic__scenario-C.png), [view](figures/figure_8_dataset-synthetic__scenario-D.png)


## Interpretation limits

The policy is not guaranteed to outperform Probability-of-Hit or other baselines. Negative and mixed outcomes remain in the saved tables. A high transcriptional-response score does not establish therapeutic usefulness, and perturbation-induced expression changes do not identify direct causal regulatory edges. STRING communities are prior association groups; control-feature clusters are artificial feature clusters. Projected feature geometry has no direct causal interpretation.

K562–RPE1 comparison requires measured outputs from both contexts under frozen method hyperparameters; missing context results are not replaced with synthetic values. Offline oracle rankings are full-information upper-bound references and are not deployable policies.


Sources and protocol are documented in the repository's `reports/bibliography.md` and `reports/research_report.md`.
