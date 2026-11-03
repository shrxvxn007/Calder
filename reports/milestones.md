# Implementation record

2026-11-03. The original workspace was empty. Local validation uses Python 3.11.15 and the exact versions in `uv.lock`.

1. Repository/configuration, four known-truth synthetic scenarios and query-only oracle: implemented. Synthetic tests verify deterministic generation, exact linear truth in scenario A, label/feature separation and duplicate query rejection.
2. Surrogates: implemented mean, ridge, fixed-prior Bayesian ridge, random Fourier Bayesian model, optional GP and five-member CPU MLP bootstrap ensemble. Analytic posterior and joint Thompson sampling checks pass.
3. Policies/calibration: all eight policies implemented under common batch/cost constraints, disjoint pilot isotonic calibration with explicit fallback and exact Gaussian-posterior Thompson sampling for Bayesian models.
4. Synthetic replay/leakage: end-to-end paired runner and pre-query snapshots implemented. Hidden-label permutation and held-out-calibration exclusion tests pass. Initial integrated suite: 68 passed. Full reference experiments and their measured outcomes are recorded separately below/in the research report.
5. Real-data ingestion: backed sparse AnnData, control-only transformations, inverse-variance weighted effects, seeded cell bootstrap and versioned network TSV supported. Toy fixtures validate the pipeline; they are not biological evidence.
6. K562 benchmark: implementation ready; actual run pending verified biological source. The public source verification attempts returned HTTP 400 for both expected objects and GCS metadata endpoints in this environment. See `data/manifests/source_verification.json`.
7. RPE1 validation: same data limitation. Frozen-hyperparameter context runner/comparison and controlled component/sensitivity studies are implemented. Synthetic ablations can run independently.
8. Saved-result figures/reports and read-only optional demo implemented. Missing biological plots are explicitly unavailable; no values are invented.

Git commits record tested working states. Wall-clock runtimes and environment details are saved with experiments; timing is excluded from deterministic trajectory equivalence.
