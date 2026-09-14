library(rjags)
library(coda)

source("01_data_priors.R")

if (!file.exists("herd_primary_summary.csv")) {
  stop("Run 04_run_models.R first")
}

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

fit_model <- function(model_file, data, inits, pars, chains=4, adapt=5000, burn=10000, iter=100000) {
  m <- jags.model(model_file, data=data, inits=inits, n.chains=chains, n.adapt=adapt)
  update(m, burn)
  coda.samples(m, pars, n.iter=iter, thin=5)
}

herd_summary <- read.csv("herd_primary_summary.csv")

get_prior <- function(parameter) {
  x <- herd_summary[herd_summary$parameter == parameter, ]
  beta_from_moments(x$mean, x$sd)
}

transfer <- list(
  se_elisa = get_prior("SeS"),
  sp_elisa = get_prior("SpS"),
  se_copro = get_prior("SeC"),
  sp_copro = get_prior("SpC")
)

widen <- function(x) {
  if (all(x == c(1,1))) x else x/5
}

make_herd_priors <- function(scenario) {
  p <- herd_priors
  tests <- c("se_elisa","sp_elisa","se_copro","sp_copro")
  prev <- c("p_arab","p_goran","p_fulani")

  if (scenario == "widened") p <- lapply(p, widen)
  if (scenario == "uniform_all") p <- lapply(p, function(x) c(1,1))
  if (scenario == "beta22_tests") p[tests] <- lapply(p[tests], function(x) c(2,2))
  if (scenario == "uniform_tests") p[tests] <- lapply(p[tests], function(x) c(1,1))
  if (scenario == "uniform_prevalence") p[prev] <- lapply(p[prev], function(x) c(1,1))
  p
}

make_slaughter_priors <- function(scenario) {
  p <- list(
    prevalence = slaughter_priors$prevalence,
    se_liver = slaughter_priors$se_liver,
    sp_liver = slaughter_priors$sp_liver,
    se_elisa = transfer$se_elisa,
    sp_elisa = transfer$sp_elisa,
    se_copro = transfer$se_copro,
    sp_copro = transfer$sp_copro
  )

  tests <- c("se_liver","sp_liver","se_elisa","sp_elisa","se_copro","sp_copro")

  if (scenario == "widened") p <- lapply(p, widen)
  if (scenario == "uniform_all") p <- lapply(p, function(x) c(1,1))
  if (scenario == "beta22_tests") p[tests] <- lapply(p[tests], function(x) c(2,2))
  if (scenario == "uniform_tests") p[tests] <- lapply(p[tests], function(x) c(1,1))
  if (scenario == "uniform_prevalence") p$prevalence <- c(1,1)
  p
}

herd_inits <- list(
  list(p=c(.08,.015,.65), SeS=.91, SpS=.94, SeC=.15, SpC=.90),
  list(p=c(.08,.015,.65), SeS=.97, SpS=.86, SeC=.40, SpC=.85),
  list(p=c(.08,.015,.65), SeS=.75, SpS=.98, SeC=.65, SpC=.95),
  list(p=c(.08,.015,.65), SeS=.85, SpS=.90, SeC=.25, SpC=.88)
)

slaughter_inits <- list(
  list(p=.08, SeL=.30, SpL=.93, SeS=.91, SpS=.94, SeC=.15, SpC=.90),
  list(p=.20, SeL=.60, SpL=.85, SeS=.97, SpS=.86, SeC=.40, SpC=.85),
  list(p=.12, SeL=.20, SpL=.97, SeS=.75, SpS=.98, SeC=.65, SpC=.95),
  list(p=.35, SeL=.45, SpL=.90, SeS=.85, SpS=.90, SeC=.25, SpC=.88)
)


zero_covariance_position <- function(a, b) {
  lower <- max(-a*b, -(1-a)*(1-b))
  upper <- min(a*(1-b), (1-a)*b)
  -lower/(upper-lower)
}

herd_dependence_inits <- lapply(herd_inits, function(z) {
  z$uSe <- zero_covariance_position(z$SeS, z$SeC)
  z$uSp <- zero_covariance_position(z$SpS, z$SpC)
  z
})

slaughter_dependence_inits <- lapply(slaughter_inits, function(z) {
  z$uSeLC <- zero_covariance_position(z$SeL, z$SeC)
  z$uSeLS <- zero_covariance_position(z$SeL, z$SeS)
  z$uSeCS <- zero_covariance_position(z$SeC, z$SeS)
  z$uSpLC <- zero_covariance_position(z$SpL, z$SpC)
  z$uSpLS <- zero_covariance_position(z$SpL, z$SpS)
  z$uSpCS <- zero_covariance_position(z$SpC, z$SpS)
  z
})

scenarios <- c("widened","uniform_all","beta22_tests","uniform_tests","uniform_prevalence")

for (s in scenarios) {
  p <- make_herd_priors(s)

  d <- list(
    y=herd_counts, n=as.integer(herd_n), oneOrient=1L,
    ap1=p$p_arab[1], bp1=p$p_arab[2],
    ap2=p$p_goran[1], bp2=p$p_goran[2],
    ap3=p$p_fulani[1], bp3=p$p_fulani[2],
    aSeS=p$se_elisa[1], bSeS=p$se_elisa[2],
    aSpS=p$sp_elisa[1], bSpS=p$sp_elisa[2],
    aSeC=p$se_copro[1], bSeC=p$se_copro[2],
    aSpC=p$sp_copro[1], bSpC=p$sp_copro[2]
  )

  x <- fit_model("02_herd_model.jags", d, herd_inits, c("p","SeS","SpS","SeC","SpC"))
  write.csv(summarise_model(x), paste0("herd_",s,".csv"), row.names=FALSE)
}

for (s in scenarios) {
  p <- make_slaughter_priors(s)

  d <- list(
    y=as.integer(slaughter_counts), n=as.integer(slaughter_n), oneOrient=1L,
    ap=p$prevalence[1], bp=p$prevalence[2],
    aSeL=p$se_liver[1], bSeL=p$se_liver[2],
    aSpL=p$sp_liver[1], bSpL=p$sp_liver[2],
    lowerSpL=0,
    aSeC=p$se_copro[1], bSeC=p$se_copro[2],
    aSpC=p$sp_copro[1], bSpC=p$sp_copro[2],
    aSeS=p$se_elisa[1], bSeS=p$se_elisa[2],
    aSpS=p$sp_elisa[1], bSpS=p$sp_elisa[2]
  )

  x <- fit_model("03_slaughterhouse_model.jags", d, slaughter_inits,
                 c("p","SeL","SpL","SeC","SpC","SeS","SpS"))
  write.csv(summarise_model(x), paste0("slaughter_",s,".csv"), row.names=FALSE)
}

p <- make_slaughter_priors("baseline")
d <- list(
  y=as.integer(slaughter_counts), n=as.integer(slaughter_n), oneOrient=1L,
  ap=p$prevalence[1], bp=p$prevalence[2],
  aSeL=p$se_liver[1], bSeL=p$se_liver[2],
  aSpL=p$sp_liver[1], bSpL=p$sp_liver[2],
  lowerSpL=.9000001,
  aSeC=p$se_copro[1], bSeC=p$se_copro[2],
  aSpC=p$sp_copro[1], bSpC=p$sp_copro[2],
  aSeS=p$se_elisa[1], bSeS=p$se_elisa[2],
  aSpS=p$sp_elisa[1], bSpS=p$sp_elisa[2]
)
x <- fit_model("03_slaughterhouse_model.jags", d, slaughter_inits,
               c("p","SeL","SpL","SeC","SpC","SeS","SpS"))
write.csv(summarise_model(x), "slaughter_liver_sp90.csv", row.names=FALSE)

# Conditional-dependence sensitivity

p <- make_herd_priors("baseline")
d <- list(
  y=herd_counts, n=as.integer(herd_n), oneOrient=1L,
  ap1=p$p_arab[1], bp1=p$p_arab[2],
  ap2=p$p_goran[1], bp2=p$p_goran[2],
  ap3=p$p_fulani[1], bp3=p$p_fulani[2],
  aSeS=p$se_elisa[1], bSeS=p$se_elisa[2],
  aSpS=p$sp_elisa[1], bSpS=p$sp_elisa[2],
  aSeC=p$se_copro[1], bSeC=p$se_copro[2],
  aSpC=p$sp_copro[1], bSpC=p$sp_copro[2]
)
x <- fit_model("06_herd_dependence.jags", d, herd_dependence_inits,
               c("p","SeS","SpS","SeC","SpC","covSe","covSp"))
write.csv(summarise_model(x), "herd_dependence.csv", row.names=FALSE)

p <- make_slaughter_priors("baseline")
d <- list(
  y=as.integer(slaughter_counts), n=as.integer(slaughter_n),
  oneOrient=1L, oneD=1L, oneN=1L,
  ap=p$prevalence[1], bp=p$prevalence[2],
  aSeL=p$se_liver[1], bSeL=p$se_liver[2],
  aSpL=p$sp_liver[1], bSpL=p$sp_liver[2],
  aSeC=p$se_copro[1], bSeC=p$se_copro[2],
  aSpC=p$sp_copro[1], bSpC=p$sp_copro[2],
  aSeS=p$se_elisa[1], bSeS=p$se_elisa[2],
  aSpS=p$sp_elisa[1], bSpS=p$sp_elisa[2]
)
x <- fit_model("07_slaughterhouse_dependence.jags", d, slaughter_dependence_inits,
               c("p","SeL","SpL","SeC","SpC","SeS","SpS",
                 "covSeLC","covSeLS","covSeCS","covSpLC","covSpLS","covSpCS"))
write.csv(summarise_model(x), "slaughter_dependence.csv", row.names=FALSE)
