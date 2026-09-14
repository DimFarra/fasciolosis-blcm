# Data and priors

# Herd: ELISA-/copro-, ELISA+/copro-, ELISA-/copro+, ELISA+/copro+
herd_counts <- rbind(
  Arab   = c(67, 2, 6, 1),
  Goran  = c(22, 0, 1, 0),
  Fulani = c(32, 78, 1, 7)
)

herd_n <- rowSums(herd_counts)

# Slaughterhouse: liver/coprology/ELISA
# ---, +--, -+-, ++-, --+, +-+, -++, +++
slaughter_counts <- c(86, 5, 4, 0, 9, 1, 3, 1)
slaughter_n <- sum(slaughter_counts)

herd_priors <- list(
  p_arab   = c(3.901, 11.037),
  p_goran  = c(1.148, 80.948),
  p_fulani = c(10.071, 11.229),
  se_elisa = c(53.581, 2.626),
  sp_elisa = c(54.706, 5.042),
  se_copro = c(4.842, 3.561),
  sp_copro = c(14.844, 4.461)
)

slaughter_priors <- list(
  prevalence = c(3.693, 10.016),
  se_liver   = c(1, 1),
  sp_liver   = c(1, 1)
)

beta_from_moments <- function(mean, sd) {
  v <- sd^2
  k <- mean * (1 - mean) / v - 1
  c(alpha = mean * k, beta = (1 - mean) * k)
}
