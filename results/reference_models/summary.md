# Exploratory development model comparison

These are actual controlled synthetic runs, not biological results. Three predeclared models were evaluated without hyperparameter tuning on 250 candidates, 24 features, four scenarios, seeds 0–2, threshold quantile 0.8, and policies random, Probability-of-Hit and CALDER. The method hyperparameters of the main benchmark remain unchanged. Three seeds are a limited development check, not final confirmatory evidence.

The random-features model uses 128 Fourier features; Bayesian ridge uses the original features. The ensemble uses five bootstrap scikit-learn MLPs with one hidden layer of 32 units and max_iter=300; disagreement is heuristic epistemic uncertainty and does not incorporate per-observation noise weights. This is a five-member MLP ensemble, not a claim to reproduce an external deep-ensemble protocol. All methods receive the same permitted candidate information and pilot split.

Means below summarize three seeds separately within each scenario/model/policy. Confirmed hits use noisy measured log effects; latent hits are a separate ground-truth diagnostic. Group coverage uses the original fixed groups. AUC is the left-step hit × cost integral beyond the pilot, including pilot hits as the initial discovery count. Every trajectory spent 50 additional cost units and selected 100 candidates including the common pilot.

| scenario | model | policy | confirmed_hits | latent_hits | group_coverage | auc_hits | brier_score | mae |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | bayesian_ridge | calder | 44.67 | 47.67 | 3 | 1330 | 0.1482 | 0.2256 |
| A | bayesian_ridge | probability_of_hit | 48 | 50.67 | 2.667 | 1516 | 0.08168 | 0.2254 |
| A | bayesian_ridge | random | 22.67 | 22.67 | 3 | 810.3 | 0.03768 | 0.2283 |
| A | ensemble | calder | 42.67 | 45.67 | 2.667 | 1329 | 0.1917 | 0.3142 |
| A | ensemble | probability_of_hit | 46 | 48.67 | 2.667 | 1475 | 0.1097 | 0.3204 |
| A | ensemble | random | 22.67 | 22.67 | 3 | 810.3 | 0.05802 | 0.3261 |
| A | random_features | calder | 32.33 | 33 | 2.667 | 958.3 | 0.312 | 0.6613 |
| A | random_features | probability_of_hit | 37.67 | 39 | 2.667 | 1228 | 0.2277 | 0.6848 |
| A | random_features | random | 22.67 | 22.67 | 3 | 810.3 | 0.1437 | 0.6578 |
| B | bayesian_ridge | calder | 18.33 | 19 | 5 | 663 | 0.2814 | 0.761 |
| B | bayesian_ridge | probability_of_hit | 27.67 | 27 | 5 | 844 | 0.3899 | 0.8113 |
| B | bayesian_ridge | random | 20 | 17 | 5 | 740.3 | 0.2 | 0.5859 |
| B | ensemble | calder | 20 | 20 | 5 | 716 | 0.232 | 0.5648 |
| B | ensemble | probability_of_hit | 26 | 25.67 | 4.333 | 853 | 0.2543 | 0.5555 |
| B | ensemble | random | 20 | 17 | 5 | 740.3 | 0.1574 | 0.5222 |
| B | random_features | calder | 21.33 | 20.67 | 5 | 748 | 0.3331 | 0.6503 |
| B | random_features | probability_of_hit | 24.67 | 22 | 4.667 | 792.3 | 0.2476 | 0.6137 |
| B | random_features | random | 20 | 17 | 5 | 740.3 | 0.1727 | 0.5267 |
| C | bayesian_ridge | calder | 21 | 18.33 | 5 | 763 | 0.2995 | 0.7073 |
| C | bayesian_ridge | probability_of_hit | 24 | 21 | 5 | 835.7 | 0.3318 | 0.6844 |
| C | bayesian_ridge | random | 24.67 | 19.33 | 5 | 860.3 | 0.2501 | 0.6694 |
| C | ensemble | calder | 21.67 | 18 | 5 | 734.7 | 0.3442 | 0.5446 |
| C | ensemble | probability_of_hit | 26 | 21 | 5 | 837.7 | 0.241 | 0.5263 |
| C | ensemble | random | 24.67 | 19.33 | 5 | 860.3 | 0.2406 | 0.5654 |
| C | random_features | calder | 24 | 21.33 | 5 | 806.7 | 0.2639 | 0.536 |
| C | random_features | probability_of_hit | 25.67 | 22.67 | 5 | 856.3 | 0.2242 | 0.455 |
| C | random_features | random | 24.67 | 19.33 | 5 | 860.3 | 0.2111 | 0.5609 |
| D | bayesian_ridge | calder | 22 | 18.33 | 5 | 684.3 | 0.4085 | 0.8243 |
| D | bayesian_ridge | probability_of_hit | 24.33 | 20.67 | 5 | 821 | 0.3915 | 0.9172 |
| D | bayesian_ridge | random | 19 | 15.33 | 5 | 694.3 | 0.2342 | 0.7958 |
| D | ensemble | calder | 23.67 | 18.33 | 5 | 780.7 | 0.3778 | 0.6768 |
| D | ensemble | probability_of_hit | 26.33 | 22.33 | 4.667 | 842.3 | 0.2886 | 0.672 |
| D | ensemble | random | 19 | 15.33 | 5 | 694.3 | 0.1851 | 0.646 |
| D | random_features | calder | 23 | 19 | 5 | 789.7 | 0.3221 | 0.7007 |
| D | random_features | probability_of_hit | 23.67 | 20.67 | 5 | 835 | 0.2742 | 0.7081 |
| D | random_features | random | 19 | 15.33 | 5 | 694.3 | 0.1969 | 0.6812 |

## Observed development comparisons

Probability-of-Hit had more mean confirmed hits than CALDER in all 12 scenario/model combinations. CALDER sometimes traded fewer hits for greater group coverage, while several scenarios had coverage already saturated at five groups. In scenario A, PoH with Bayesian ridge found 48.0 mean confirmed hits, the ensemble 46.0 and RFF 37.67; the model-family ranking therefore differed from simply assuming that the nonlinear surrogate must be best. These small comparisons are exploratory and do not change the predeclared main model or its settings.

| scenario | model | hit_difference_calder_minus_poh | group_difference_calder_minus_poh |
| --- | --- | --- | --- |
| A | bayesian_ridge | -3.333 | 0.3333 |
| A | ensemble | -3.333 | 0 |
| A | random_features | -5.333 | 0 |
| B | bayesian_ridge | -9.333 | 0 |
| B | ensemble | -6 | 0.6667 |
| B | random_features | -3.333 | 0.3333 |
| C | bayesian_ridge | -3 | 0 |
| C | ensemble | -4.333 | 0 |
| C | random_features | -1.667 | 0 |
| D | bayesian_ridge | -2.333 | 0 |
| D | ensemble | -2.667 | 0.3333 |
| D | random_features | -0.6667 | 0 |

## Warning record

| model | category | message | count |
| --- | --- | --- | --- |
| ensemble | ConvergenceWarning | Stochastic Optimizer: Maximum iterations (300) reached and the optimization hasn't converged yet. | 153 |

Convergence warnings are retained rather than suppressed or resolved by tuning these runs. A warning means the MLP optimizer exhausted its iteration limit; ensemble comparisons therefore reflect this fixed training budget. The full warning sequence is retained in each model directory as warnings.json and warnings.txt. No claim of ensemble convergence or model-family superiority follows from these small runs.


## Reproduction and audit

Run `.venv/bin/python results/reference_models/run_comparison.py` from the repository root after the frozen installation. The script refuses existing incomplete runs and reuses completed model runs. To create independent replicas, use the saved per-model config.json snapshots with `calder benchmark --config <snapshot> --run-dir <new-empty-directory>`. Each model has its own configuration hash, Git metadata, environment information, all seed-level outcomes, pre-query prediction snapshots, audit logs and figures. Aggregate files here concatenate those real outputs.

Repeated seeds describe algorithmic variation on one fixed synthetic pool; they are not biological replications. Calibration diagnostics on actively selected examples are selection biased. Negative/mixed results are retained; results must not be used to retune the final RPE1 benchmark.
