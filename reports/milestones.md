# Implementation record

2026-11-03. The original workspace was empty. Local validation uses Python 3.11.15 and the exact versions in `uv.lock`.

1. Repository/configuration, four known-truth synthetic scenarios and query-only oracle: implemented. Synthetic tests verify deterministic generation, exact linear truth in scenario A, label/feature separation and duplicate query rejection.
2. Surrogates: implemented mean, ridge, fixed-prior Bayesian ridge, random Fourier Bayesian model, optional GP and five-member CPU MLP bootstrap ensemble. Analytic posterior and joint Thompson sampling checks pass.
3. Policies/calibration: all eight policies implemented under common batch/cost constraints, disjoint pilot isotonic calibration with explicit fallback and exact Gaussian-posterior Thompson sampling for Bayesian models.
4. Synthetic replay/leakage: end-to-end paired runner and pre-query snapshots implemented. Hidden-label permutation and held-out-calibration exclusion tests pass. Initial integrated suite: 68 passed. Full reference experiments and their measured outcomes are recorded separately below/in the research report.
5. Real-data ingestion: backed sparse AnnData, control-only transformations, inverse-variance weighted effects, seeded cell bootstrap and versioned network TSV supported. Toy fixtures validate the pipeline; they are not biological evidence.
6. K562 benchmark: completed 270 trajectories (10 seeds × three thresholds × eight policies plus reference), using 1,855 eligible perturbations from checksum-verified scPerturb Zenodo essential-panel counts. Run: `results/reference_k562`; execution commit `7e95d25`, tracked source clean. Original scPertEval GCS access failures remain separately recorded.
7. RPE1 validation: completed 270 trajectories over 2,045 eligible perturbations with identical frozen model/acquisition/protocol and execution source hash. Run: `results/reference_rpe1`. Also completed 840 synthetic ablation, 108 development model-comparison and 216 sensitivity trajectories. No RPE1 retuning occurred.
8. Saved-result figures/reports and read-only optional demo implemented. All eight figure families have genuine supporting runs; combined biological comparison is `results/reference_contexts`. The main report preserves lower hit discovery for CALDER and small artificial-cluster coverage gains. Network-feature and alternative-effect benchmark variants are supported but have not been empirically run.

Git commits record tested working states. Wall-clock runtimes and environment details are saved with experiments; timing is excluded from deterministic trajectory equivalence.

Final verification: 78 tests passed; Ruff check and format passed; wheel/sdist build passed; a wheel installed into a separate fresh environment. The complete experiment set contains 2,784 replay trajectories, including 540 biological-context trajectories. Raw H5ADs and full prediction snapshots stay outside Git; saved tabular outcomes, summaries, figures and provenance can be replayed from configurations.
