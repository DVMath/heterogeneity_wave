# heterogeneity_wave — core pipeline

## Inputs (data files, already in repo root)

* `dataset_PNAS_SARSCOV2_655_anon_clean.csv` — main cohort (kinetics,
survival, transmission).
* `elife-69302-data1-v2-converted.csv` — eLife donor/contact cohort
(kinetics + transmission-risk likelihood).
* `FR_IDF.csv` — Île-de-France population case/death counts (SIR).

### Step 1 — Fit the joint model

```r
Rscript run_elife_expdecay_concurrent_model.R
```

* Model file: `myclaphammodel_v17_elife_expdecay_concurrent.stan`.
* Loads the three CSVs above, builds `mydata`, fits with 4 chains,
`adapt_delta = 0.99`, `max_treedepth = 15`, warmup 4000 / iter 5000.
* Writes, with tag `expdecay_concurrent_lag0`:

  * `elife_kinetics_fit_expdecay_concurrent_lag0.RData` (full `stanfit`)
  * `elife_kinetics_posterior_expdecay_concurrent_lag0.RData` (extracted
posterior draws — this is what Step 2 actually loads)
  * `elife_kinetics_summary_expdecay_concurrent_lag0.RData`
  * a handful of diagnostic `mcmc_pairs` PNGs

### Step 2 — Post-processing: all tables and figures

```r
Rscript post_analysis_v17_expdecay_concurrent.R
```

* Loads `elife_kinetics_posterior_expdecay_concurrent_lag0.RData` from
Step 1.
* Regenerates Table 1 (grouped Kinetics&Survival / Transmission /
Population), Table 2 (grouped Kinetics / Survival / Transmission /
Population, plain-language labels), Table 3, Table 4, and Figures
2 / 3 / S1 / S2, all written to `figures_new/`.
* Runs two Monte-Carlo simulation loops (the slow part) before Figure 3,
then **saves a checkpoint**: `checkpoint_before_fig3.RData`. Everything after that point — Figure 3's
plotting code — is what Step 2b re-runs quickly.

### Step 2b — Fast Figure 3 iteration

Only needed to tweak Figure 3 (colors, fonts, labels)
without re-running Step 2 from scratch:

```r
Rscript post_analysis_v17_expdecay_figure3_only.R
```

* Requires `checkpoint_before_fig3.RData` from Step 2 to already exist.
* Loads it and re-executes only the Figure 3 plotting block, writing the
updated `figures_new/sir_extreme_individuals.pdf/.tiff`.
