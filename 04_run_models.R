library(rjags)
library(coda)

source("01_data_priors.R")

summarise_model <- function(x) {
  x <- as.matrix(x)
  q <- t(apply(x, 2, quantile, c(.025, .5, .975)))

  data.frame(
    parameter = colnames(x),
    lower = q[,1],
    median = q[,2],
    upper = q[,3],
    mean = colMeans(x),
    sd = apply(x, 2, sd),
    row.names = NULL
  )
}

# Herd model

herd_data <- list(
  y = herd_counts,
  n = as.integer(herd_n),
  oneOrient = 1L,
  ap1 = herd_priors$p_arab[1],
  bp1 = herd_priors$p_arab[2],
  ap2 = herd_priors$p_goran[1],
  bp2 = herd_priors$p_goran[2],
  ap3 = herd_priors$p_fulani[1],
  bp3 = herd_priors$p_fulani[2],
  aSeS = herd_priors$se_elisa[1],
  bSeS = herd_priors$se_elisa[2],
  aSpS = herd_priors$sp_elisa[1],
  bSpS = herd_priors$sp_elisa[2],
  aSeC = herd_priors$se_copro[1],
  bSeC = herd_priors$se_copro[2],
  aSpC = herd_priors$sp_copro[1],
  bSpC = herd_priors$sp_copro[2]
)

herd_inits <- list(
  list(p=c(.08,.015,.65), SeS=.91, SpS=.94, SeC=.15, SpC=.90,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=910261),
  list(p=c(.08,.015,.65), SeS=.97, SpS=.86, SeC=.40, SpC=.85,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=910262),
  list(p=c(.08,.015,.65), SeS=.75, SpS=.98, SeC=.65, SpC=.95,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=910263),
  list(p=c(.08,.015,.65), SeS=.85, SpS=.90, SeC=.25, SpC=.88,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=910264)
)

herd_model <- jags.model(
  "02_herd_model.jags",
  data = herd_data,
  inits = herd_inits,
  n.chains = 4,
  n.adapt = 5000
)

update(herd_model, 10000)

herd_draws <- coda.samples(
  herd_model,
  c("p","SeS","SpS","SeC","SpC"),
  n.iter = 100000,
  thin = 5
)

herd_summary <- summarise_model(herd_draws)
write.csv(herd_summary, "herd_primary_summary.csv", row.names = FALSE)


# Priors for slaughterhouse model

get_prior <- function(parameter) {
  x <- herd_summary[herd_summary$parameter == parameter, ]
  beta_from_moments(x$mean, x$sd)
}

transfer <- rbind(
  SeS = get_prior("SeS"),
  SpS = get_prior("SpS"),
  SeC = get_prior("SeC"),
  SpC = get_prior("SpC")
)

write.csv(transfer, "transferred_test_priors.csv")


# Slaughterhouse model

slaughter_data <- list(
  y = as.integer(slaughter_counts),
  n = as.integer(slaughter_n),
  oneOrient = 1L,
  ap = slaughter_priors$prevalence[1],
  bp = slaughter_priors$prevalence[2],
  aSeL = slaughter_priors$se_liver[1],
  bSeL = slaughter_priors$se_liver[2],
  aSpL = slaughter_priors$sp_liver[1],
  bSpL = slaughter_priors$sp_liver[2],\n  lowerSpL = 0,
  aSeC = transfer["SeC","alpha"],
  bSeC = transfer["SeC","beta"],
  aSpC = transfer["SpC","alpha"],
  bSpC = transfer["SpC","beta"],
  aSeS = transfer["SeS","alpha"],
  bSeS = transfer["SeS","beta"],
  aSpS = transfer["SpS","alpha"],
  bSpS = transfer["SpS","beta"]
)

slaughter_inits <- list(
  list(p=.08, SeL=.30, SpL=.93, SeS=.91, SpS=.94, SeC=.15, SpC=.90,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=912261),
  list(p=.20, SeL=.60, SpL=.85, SeS=.97, SpS=.86, SeC=.40, SpC=.85,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=912262),
  list(p=.12, SeL=.20, SpL=.97, SeS=.75, SpS=.98, SeC=.65, SpC=.95,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=912263),
  list(p=.35, SeL=.45, SpL=.90, SeS=.85, SpS=.90, SeC=.25, SpC=.88,
       .RNG.name="base::Mersenne-Twister", .RNG.seed=912264)
)

slaughter_model <- jags.model(
  "03_slaughterhouse_model.jags",
  data = slaughter_data,
  inits = slaughter_inits,
  n.chains = 4,
  n.adapt = 5000
)

update(slaughter_model, 10000)

slaughter_draws <- coda.samples(
  slaughter_model,
  c("p","SeL","SpL","SeC","SpC","SeS","SpS"),
  n.iter = 100000,
  thin = 5
)

slaughter_summary <- summarise_model(slaughter_draws)
write.csv(slaughter_summary, "slaughterhouse_primary_summary.csv", row.names = FALSE)
