# Exploratory one-factor sensitivity experiments

These are 216 actual controlled synthetic trajectories from nine predeclared development variants. The frozen RFF model, features and acquisition weights are unchanged. Four scenarios, seeds 0–2, threshold quantile 0.8 and the policies Probability-of-Hit and CALDER were used. No main or biological hyperparameter was selected from these outcomes.

The one-factor grid varies pilot fraction (0.1, 0.2, 0.3), maximum batch size (5, 10, 20), and simulated observation-noise scale (0.05, 0.2, 0.5). All other fields retain the supplied base config. The base setting appears once in each factor axis and is intentionally retained, not counted as an independent population. Its three independent reruns produced exactly identical scientific selections and observed outcomes, checked automatically.

Compare policies within the same variant and scenario. Pilot-size changes also change total measured candidates and the pilot-derived threshold. Noise changes affect noisy labels and their pilot-derived threshold. Consequently, unadjusted hit counts across variants do not measure one common fixed-threshold estimand. Confirmed hits include the pilot; active-only and latent-hit counts are saved separately.

For final measured hits, CALDER−PoH mean differences were negative in 34, zero in 0 and positive in 2 of the 36 scenario/variant cells. These descriptive cell counts are not independent statistical replications. No cross-scenario confidence interval is computed.

## Means by scenario and variant

| variant | scenario | policy | confirmed_hits | active_hits | group_coverage | brier_score | calibration_fitted_fraction | n_selected |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| batch_size_10 | A | calder | 32.33 | 20.33 | 2.667 | 0.312 | 1 | 100 |
| batch_size_10 | A | probability_of_hit | 37.67 | 25.67 | 2.667 | 0.2277 | 1 | 100 |
| batch_size_10 | B | calder | 21.33 | 11 | 5 | 0.3331 | 1 | 100 |
| batch_size_10 | B | probability_of_hit | 24.67 | 14.33 | 4.667 | 0.2476 | 1 | 100 |
| batch_size_10 | C | calder | 24 | 13 | 5 | 0.2639 | 1 | 100 |
| batch_size_10 | C | probability_of_hit | 25.67 | 14.67 | 5 | 0.2242 | 1 | 100 |
| batch_size_10 | D | calder | 23 | 12.67 | 5 | 0.3221 | 1 | 100 |
| batch_size_10 | D | probability_of_hit | 23.67 | 13.33 | 5 | 0.2742 | 1 | 100 |
| batch_size_20 | A | calder | 29 | 17 | 2.667 | 0.3335 | 1 | 100 |
| batch_size_20 | A | probability_of_hit | 36.33 | 24.33 | 2.667 | 0.2218 | 1 | 100 |
| batch_size_20 | B | calder | 20 | 9.667 | 4.667 | 0.3024 | 1 | 100 |
| batch_size_20 | B | probability_of_hit | 25 | 14.67 | 4.333 | 0.2506 | 1 | 100 |
| batch_size_20 | C | calder | 22.67 | 11.67 | 5 | 0.2595 | 1 | 100 |
| batch_size_20 | C | probability_of_hit | 26.67 | 15.67 | 5 | 0.2259 | 1 | 100 |
| batch_size_20 | D | calder | 22 | 11.67 | 5 | 0.2717 | 1 | 100 |
| batch_size_20 | D | probability_of_hit | 23 | 12.67 | 5 | 0.27 | 1 | 100 |
| batch_size_5 | A | calder | 30.67 | 18.67 | 2.667 | 0.3387 | 1 | 100 |
| batch_size_5 | A | probability_of_hit | 37.33 | 25.33 | 2.667 | 0.2213 | 1 | 100 |
| batch_size_5 | B | calder | 19.67 | 9.333 | 5 | 0.2945 | 1 | 100 |
| batch_size_5 | B | probability_of_hit | 23.67 | 13.33 | 4.667 | 0.2422 | 1 | 100 |
| batch_size_5 | C | calder | 24 | 13 | 5 | 0.3262 | 1 | 100 |
| batch_size_5 | C | probability_of_hit | 26.33 | 15.33 | 5 | 0.2262 | 1 | 100 |
| batch_size_5 | D | calder | 20.33 | 10 | 4.667 | 0.3145 | 1 | 100 |
| batch_size_5 | D | probability_of_hit | 22.33 | 12 | 4.667 | 0.2587 | 1 | 100 |
| noise_scale_0.05 | A | calder | 32.33 | 20.33 | 3 | 0.3449 | 1 | 100 |
| noise_scale_0.05 | A | probability_of_hit | 39.67 | 27.67 | 3 | 0.2355 | 1 | 100 |
| noise_scale_0.05 | B | calder | 24 | 12.33 | 4.333 | 0.3761 | 1 | 100 |
| noise_scale_0.05 | B | probability_of_hit | 27.33 | 15.67 | 4.667 | 0.2761 | 1 | 100 |
| noise_scale_0.05 | C | calder | 24 | 13.33 | 5 | 0.3547 | 1 | 100 |
| noise_scale_0.05 | C | probability_of_hit | 26.67 | 16 | 4.667 | 0.2519 | 1 | 100 |
| noise_scale_0.05 | D | calder | 23.33 | 13 | 5 | 0.264 | 1 | 100 |
| noise_scale_0.05 | D | probability_of_hit | 26.33 | 16 | 5 | 0.2983 | 1 | 100 |
| noise_scale_0.2 | A | calder | 32.33 | 20.33 | 2.667 | 0.312 | 1 | 100 |
| noise_scale_0.2 | A | probability_of_hit | 37.67 | 25.67 | 2.667 | 0.2277 | 1 | 100 |
| noise_scale_0.2 | B | calder | 21.33 | 11 | 5 | 0.3331 | 1 | 100 |
| noise_scale_0.2 | B | probability_of_hit | 24.67 | 14.33 | 4.667 | 0.2476 | 1 | 100 |
| noise_scale_0.2 | C | calder | 24 | 13 | 5 | 0.2639 | 1 | 100 |
| noise_scale_0.2 | C | probability_of_hit | 25.67 | 14.67 | 5 | 0.2242 | 1 | 100 |
| noise_scale_0.2 | D | calder | 23 | 12.67 | 5 | 0.3221 | 1 | 100 |
| noise_scale_0.2 | D | probability_of_hit | 23.67 | 13.33 | 5 | 0.2742 | 1 | 100 |
| noise_scale_0.5 | A | calder | 31.67 | 19.33 | 3.333 | 0.2769 | 1 | 100 |
| noise_scale_0.5 | A | probability_of_hit | 37.67 | 25.33 | 3 | 0.2073 | 1 | 100 |
| noise_scale_0.5 | B | calder | 22.67 | 12.33 | 4.667 | 0.2666 | 1 | 100 |
| noise_scale_0.5 | B | probability_of_hit | 26 | 15.67 | 4.333 | 0.2378 | 1 | 100 |
| noise_scale_0.5 | C | calder | 23 | 12.67 | 5 | 0.3583 | 1 | 100 |
| noise_scale_0.5 | C | probability_of_hit | 24 | 13.67 | 5 | 0.2028 | 1 | 100 |
| noise_scale_0.5 | D | calder | 22.33 | 11.33 | 5 | 0.2819 | 0.6667 | 100 |
| noise_scale_0.5 | D | probability_of_hit | 21.33 | 10.33 | 4.667 | 0.187 | 0.6667 | 100 |
| pilot_fraction_0.1 | A | calder | 20.33 | 15 | 2.333 | 0.13 | 0 | 75 |
| pilot_fraction_0.1 | A | probability_of_hit | 30.67 | 25.33 | 2.333 | 0.2236 | 0 | 75 |
| pilot_fraction_0.1 | B | calder | 18.33 | 13 | 4.333 | 0.1936 | 0 | 75 |
| pilot_fraction_0.1 | B | probability_of_hit | 18.67 | 13.33 | 3.667 | 0.2107 | 0 | 75 |
| pilot_fraction_0.1 | C | calder | 20 | 14.33 | 5 | 0.1968 | 0 | 75 |
| pilot_fraction_0.1 | C | probability_of_hit | 19.33 | 13.67 | 4.667 | 0.2043 | 0 | 75 |
| pilot_fraction_0.1 | D | calder | 17.33 | 12.33 | 5 | 0.1932 | 0 | 75 |
| pilot_fraction_0.1 | D | probability_of_hit | 20 | 15 | 5 | 0.2165 | 0 | 75 |
| pilot_fraction_0.2 | A | calder | 32.33 | 20.33 | 2.667 | 0.312 | 1 | 100 |
| pilot_fraction_0.2 | A | probability_of_hit | 37.67 | 25.67 | 2.667 | 0.2277 | 1 | 100 |
| pilot_fraction_0.2 | B | calder | 21.33 | 11 | 5 | 0.3331 | 1 | 100 |
| pilot_fraction_0.2 | B | probability_of_hit | 24.67 | 14.33 | 4.667 | 0.2476 | 1 | 100 |
| pilot_fraction_0.2 | C | calder | 24 | 13 | 5 | 0.2639 | 1 | 100 |
| pilot_fraction_0.2 | C | probability_of_hit | 25.67 | 14.67 | 5 | 0.2242 | 1 | 100 |
| pilot_fraction_0.2 | D | calder | 23 | 12.67 | 5 | 0.3221 | 1 | 100 |
| pilot_fraction_0.2 | D | probability_of_hit | 23.67 | 13.33 | 5 | 0.2742 | 1 | 100 |
| pilot_fraction_0.3 | A | calder | 30.67 | 14.33 | 2.667 | 0.3137 | 1 | 125 |
| pilot_fraction_0.3 | A | probability_of_hit | 36.33 | 20 | 2.333 | 0.2113 | 1 | 125 |
| pilot_fraction_0.3 | B | calder | 37.67 | 18.67 | 5 | 0.3621 | 1 | 125 |
| pilot_fraction_0.3 | B | probability_of_hit | 39.67 | 20.67 | 5 | 0.2713 | 1 | 125 |
| pilot_fraction_0.3 | C | calder | 30.33 | 15 | 5 | 0.287 | 1 | 125 |
| pilot_fraction_0.3 | C | probability_of_hit | 31.33 | 16 | 5 | 0.2485 | 1 | 125 |
| pilot_fraction_0.3 | D | calder | 25.33 | 10.67 | 5 | 0.2814 | 1 | 125 |
| pilot_fraction_0.3 | D | probability_of_hit | 27.33 | 12.67 | 5 | 0.2623 | 1 | 125 |

## Paired measured-hit differences

Positive values favour CALDER. Intervals bootstrap whole paired seeds separately in each scenario and variant; three seeds give a limited exploratory interval, conditional on this fixed synthetic pool. Other metric contrasts are in paired_differences.csv. Multiple contrasts are descriptive, with no adjusted superiority claim.

| variant | scenario | mean_difference | ci_lower | ci_upper | n_pairs |
| --- | --- | --- | --- | --- | --- |
| batch_size_10 | A | -5.333 | -9 | -3 | 3 |
| batch_size_10 | B | -3.333 | -4 | -2 | 3 |
| batch_size_10 | C | -1.667 | -3 | -1 | 3 |
| batch_size_10 | D | -0.6667 | -2 | 1 | 3 |
| batch_size_20 | A | -7.333 | -13 | -4 | 3 |
| batch_size_20 | B | -5 | -9 | 0 | 3 |
| batch_size_20 | C | -4 | -5 | -3 | 3 |
| batch_size_20 | D | -1 | -1 | -1 | 3 |
| batch_size_5 | A | -6.667 | -11 | -2 | 3 |
| batch_size_5 | B | -4 | -6 | -2 | 3 |
| batch_size_5 | C | -2.333 | -4 | -1 | 3 |
| batch_size_5 | D | -2 | -3 | -1 | 3 |
| noise_scale_0.05 | A | -7.333 | -10 | -4 | 3 |
| noise_scale_0.05 | B | -3.333 | -6 | 0 | 3 |
| noise_scale_0.05 | C | -2.667 | -3 | -2 | 3 |
| noise_scale_0.05 | D | -3 | -7 | -1 | 3 |
| noise_scale_0.2 | A | -5.333 | -9 | -3 | 3 |
| noise_scale_0.2 | B | -3.333 | -4 | -2 | 3 |
| noise_scale_0.2 | C | -1.667 | -3 | -1 | 3 |
| noise_scale_0.2 | D | -0.6667 | -2 | 1 | 3 |
| noise_scale_0.5 | A | -6 | -9 | -4 | 3 |
| noise_scale_0.5 | B | -3.333 | -5 | -1 | 3 |
| noise_scale_0.5 | C | -1 | -6 | 2 | 3 |
| noise_scale_0.5 | D | 1 | -2 | 4 | 3 |
| pilot_fraction_0.1 | A | -10.33 | -18 | -5 | 3 |
| pilot_fraction_0.1 | B | -0.3333 | -5 | 2 | 3 |
| pilot_fraction_0.1 | C | 0.6667 | -2 | 3 | 3 |
| pilot_fraction_0.1 | D | -2.667 | -4 | -1 | 3 |
| pilot_fraction_0.2 | A | -5.333 | -9 | -3 | 3 |
| pilot_fraction_0.2 | B | -3.333 | -4 | -2 | 3 |
| pilot_fraction_0.2 | C | -1.667 | -3 | -1 | 3 |
| pilot_fraction_0.2 | D | -0.6667 | -2 | 1 | 3 |
| pilot_fraction_0.3 | A | -5.667 | -8 | -4 | 3 |
| pilot_fraction_0.3 | B | -2 | -2 | -2 | 3 |
| pilot_fraction_0.3 | C | -1 | -5 | 3 | 3 |
| pilot_fraction_0.3 | D | -2 | -6 | 2 | 3 |

## Warnings and limitations

| variant | category | message | count |
| --- | --- | --- | --- |
| pilot_fraction_0.1 | RuntimeWarning | Using raw probabilities: too few calibration examples or insufficient class/probability variation | 192 |
| noise_scale_0.5 | RuntimeWarning | Using raw probabilities: too few calibration examples or insufficient class/probability variation | 16 |

Small calibration reserves fall back to raw probabilities and retain explicit warnings. Calibration is fitted against measured pilot labels, whereas raw epistemic probabilities concern latent effects; prequential acquired-label diagnostics remain selection biased. Repeated seeds are algorithmic runs, not biological replications.


## Reproduction

From a frozen installation at repository root, run `uv run --frozen --extra dev python scripts/run_sensitivity.py --config configs/sensitivity.yaml --output-dir results/reference_sensitivity`. Existing completed variants are reused only when their exact config hashes match. Use another output directory for fresh replay. Each variant retains config.json, source hashes, Git metadata, full seed endpoints, selections, predictions and audit logs. make_report=False avoids redundant per-variant figure generation; use `calder report --run-dir <variant-directory>` when needed. No absent figure is replaced by invented values.
