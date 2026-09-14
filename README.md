# Fasciolosis BLCM

Minimal code for the primary models and sensitivity analyses reported in the manuscript.

Files:
- `01_data_priors.R`: counts and priors
- `02_herd_model.jags`: primary herd model
- `03_slaughterhouse_model.jags`: primary slaughterhouse model
- `04_run_models.R`: runs the primary models
- `05_sensitivity_analysis.R`: prior and structural sensitivity analyses
- `06_herd_dependence.jags`: herd conditional-dependence sensitivity
- `07_slaughterhouse_dependence.jags`: slaughterhouse conditional-dependence sensitivity

Requirements: R, JAGS, `rjags`, `coda`.

Run:

```r
Rscript 04_run_models.R
Rscript 05_sensitivity_analysis.R
```

The herd model is run first. Its posterior estimates for ELISA and coprology sensitivity and specificity are used to define the corresponding slaughterhouse priors.

The main models use conditional independence. Conditional dependence is included only as a separate sensitivity analysis.
