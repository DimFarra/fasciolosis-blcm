# Fasciolosis BLCM

Minimal code to reproduce the main models.

Files:
- `01_data_priors.R`: counts and priors
- `02_herd_model.jags`: herd model
- `03_slaughterhouse_model.jags`: slaughterhouse model
- `04_run_models.R`: runs the models

Requirements: R, JAGS, `rjags`, `coda`.

Run:

```r
Rscript 04_run_models.R
```

The herd model is run first. Its posterior estimates for ELISA and coprology sensitivity and specificity are used to define the corresponding slaughterhouse priors.

The shared code contains the main conditional-independence models only.
