library(tidyverse)
library(rstan)
library(bayesplot)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

# ---------------------------------------------------------------------------
# Concurrent-Vt variant of run_elife_expdecay_model.R: fits
# myclaphammodel_v17_elife_expdecay_concurrent.stan (main cohort's own
# lptrans driven by vt[jj, ll] instead of vt[jj, ll-1]) with
# trans_lag_days = 0, so both halves of the model are concurrent -- the
# exp-decay-curve counterpart to run_elife_kinetics_concurrent_model.R, run
# with identical settings so the two eLife curve families are comparable
# under the same main-cohort transmission timing convention. Motivated by
# a finding from the up-down concurrent run: for the ~small number of
# contacts exposed well before the donor's symptom onset (t_exposure as
# negative as -21), the up-down curve's poorly-identified rise-rate
# (alpha_elife, essentially unconstrained for most cases -- see
# diagnostics_elife_upslope_check.R) extrapolates to physically
# meaningless donor viral loads (as low as -37 log10 copies/mL, ~38 orders
# of magnitude below the assay's own LOQ). The decline-only exp-decay
# curve has no rise-rate parameter and so cannot produce this pathology.
# Output filenames use an "expdecay_concurrent_lag0" tag so they never
# collide with the base expdecay run's (expdecay_lag0/lag1) files.
#
# Part 1 duplicated from analysis_comb_model.R (not sourced), same as
# run_elife_kinetics_model.R.
# ---------------------------------------------------------------------------
dfi <- read_csv2("./dataset_PNAS_SARSCOV2_655_anon_clean.csv")

dfi %>% filter(!is.na(sex)) -> dfi_filt

dfi_filt %>%
  filter(type == 1) %>%
  group_by(ID) %>%
  arrange(ID, time_monolix) %>%
  mutate(y_lag = lag(y)) %>%
  mutate(diff_y = y - y_lag) %>%
  mutate(miny = min(diff_y, na.rm = TRUE)) %>%
  mutate(maxy = max(diff_y, na.rm = TRUE)) %>%
  mutate(flag = if_else(miny == 0 & maxy == 0, 1, 0)) %>%
  filter(flag == 1) -> dobsequal

dcount <- dfi_filt %>% filter(type == 1) %>% group_by(ID) %>% summarise(count = n())
xx <- dcount %>% filter(count > 3)

dfi_filt_a <- dfi_filt %>% filter(ID %in% xx$ID)
dequal <- dfi_filt_a %>% filter(ID %in% dobsequal$ID)
IDs_equal <- unique(dequal$ID)
dfi_filt_b <- dfi_filt_a %>% filter(!(ID %in% IDs_equal))
dfi_filt <- dfi_filt_b

dfi_filt %>%
  filter(type == 2) %>%
  group_by(ID) %>%
  mutate(died = max(y)) %>%
  filter(time_monolix != 0) %>%
  dplyr::select(ID, died, time_monolix) %>%
  distinct(ID, died, time_monolix) %>%
  mutate(tevent = as.integer(time_monolix)) %>%
  select(!time_monolix) -> df_fatal

dfi_filt2 <- left_join(dfi_filt, df_fatal, by = c("ID" = "ID"))

dfi_filt2 %>%
  filter(type == 1) %>%
  ungroup() %>%
  mutate(vir = 10^y) %>%
  mutate(ID_fct = as.factor(ID)) %>%
  mutate(i_pt = as.numeric(ID_fct)) %>%
  mutate(sero = 1) %>%
  mutate(is.over65 = age_cat_reg) %>%
  mutate(is.male = 1 - sex) %>%
  mutate(severe = died) %>%
  group_by(i_pt) %>%
  mutate(Day = time_monolix) %>%
  arrange(Day) %>%
  mutate(Day_lag = lag(Day)) %>%
  mutate(diff_day = if_else(is.na(Day_lag), 1, Day - Day_lag)) %>%
  mutate(prod_vir_dt = vir * diff_day) %>%
  mutate(maxi = max(vir)) %>%
  mutate(maxv = sum(prod_vir_dt, na.rm = TRUE)) %>%
  mutate(peak_vir = max(vir, na.rm = TRUE)) %>%
  mutate(max_Day = max(Day)) %>%
  mutate(min_Day = min(Day)) %>%
  mutate(days_obs = max_Day - min_Day + 1) %>%
  mutate(rate_vir = maxv / days_obs) -> dfc3i_filt

dfc3i_filt %>% filter(type == 1) -> dfc3_filt

dfc3_filt %>%
  dplyr::select(i_pt, is.male, is.over65, severe, sero, rate_vir,
                tevent, peak_vir, maxv) %>%
  distinct(i_pt, is.male, is.over65, severe, sero, rate_vir,
           tevent, peak_vir, maxv) -> dfu

# population (Ile-de-France) data
YDF <- read_csv("FR_IDF.csv")
YDFcases <- YDF[70:132, "new_confirmed"]
YDFdeaths <- YDF[70:132, "new_deceased"]
dfYDF <- data.frame(YDFcases = unlist(YDFcases), YDFdeaths = unlist(YDFdeaths),
                     tc_obs = 1:length(unlist(YDFcases)))
dfYDF %>% replace_na(list(YDFdeaths = 0)) -> dfYDF

ddiv <- 7
dfYDF %>%
  mutate(sem = 1 + (tc_obs - 1) %/% ddiv) %>%
  add_count(sem, wt = YDFcases, name = "cases7") %>%
  add_count(sem, wt = YDFdeaths, name = "deaths7") %>%
  dplyr::select(sem, cases7, deaths7) %>%
  distinct(sem, cases7, deaths7) %>%
  mutate(tc_obs = sem * ddiv) -> xxYDF

ncases <- unlist(xxYDF$cases7)
nd <- unlist(xxYDF$deaths7)
tc_obs <- unlist(xxYDF$tc_obs)
npop <- unlist(YDF[1, "population"])

mods <- lm(log(cases7) ~ tc_obs, data = xxYDF[1:3, ])
Y0rough_estim <- exp(coef(mods)[1])

gamma_fixed <- 1 / 6.5
CFRrough <- sum(nd) / sum(ncases)
mu <- CFRrough * gamma_fixed / (1 - CFRrough)
betarough <- unname(coef(mods)[2] + mu + gamma_fixed)

n_t <- max(tc_obs)
ts <- seq(1, max(tc_obs), by = 1)
t0 <- 0
n_tcobs <- length(tc_obs)

mydata <- list(
  y_obs = log(dfc3_filt$vir),
  ysev = dfu$severe,
  N_ind = max(dfc3_filt$i_pt),
  N_obs = dim(dfc3_filt)[1],
  male = dfu$is.male,
  over65 = dfu$is.over65,
  tevent = dfu$tevent,
  nserotype = 1,
  serotype = dfu$sero,
  ratevir = dfu$rate_vir,
  peak = dfu$peak_vir,
  area = dfu$maxv,
  eps0 = -8.91,
  eps1 = 1.368,
  myinf = 30,
  sigma_meas = 0.5,
  d_fixed = 0.0,
  t = dfc3_filt$Day,
  ind = dfc3_filt$i_pt,
  ts = ts,
  t0 = t0,
  ncases = ncases,
  ndeaths = nd,
  CFRrough = CFRrough,
  betarough = betarough,
  n_t = n_t,
  tc_obs = tc_obs,
  popnorm = npop,
  npop = npop,
  Y0 = Y0rough_estim,
  n_tcobs = n_tcobs
)
mydata$nresp <- max(mydata$tevent)

# ---------------------------------------------------------------------------
# Part 2: eLife contact-tracing kinetics data, identical to
# run_elife_kinetics_model.R (male_elife/over65_elife/d_fixed_elife are
# harmless unused extras in mydata for this variant -- the exp-decay Stan
# file's data block doesn't declare them).
# ---------------------------------------------------------------------------
df_elife <- read_csv("./elife-69302-data1-v2-converted.csv")

cases_elife <- df_elife %>%
  distinct(case_id, case_sex, case_age, case_earliest_date_symptoms,
           case_day0_VL, case_day0_swab_date,
           case_day3_VL, case_day3_swab_date,
           case_day7_VL, case_day7_swab_date) %>%
  mutate(
    symp_date = as.Date(case_earliest_date_symptoms),
    is.male = if_else(case_sex == 0, 1, 0),   # case_sex==1 is coded "woman" in this dataset
    is.over65 = if_else(case_age >= 65, 1, 0)
  ) %>%
  arrange(case_id) %>%
  mutate(case_idx = row_number())

N_ind_elife <- nrow(cases_elife)

vl_long <- cases_elife %>%
  transmute(case_idx,
            d0_vl = case_day0_VL, d0_t = as.numeric(as.Date(case_day0_swab_date) - symp_date),
            d3_vl = case_day3_VL, d3_t = as.numeric(as.Date(case_day3_swab_date) - symp_date),
            d7_vl = case_day7_VL, d7_t = as.numeric(as.Date(case_day7_swab_date) - symp_date)) %>%
  pivot_longer(cols = -case_idx, names_to = c("visit", ".value"), names_sep = "_") %>%
  filter(!is.na(vl), !is.na(t)) %>%
  arrange(case_idx, t) %>%
  group_by(case_idx) %>%
  mutate(t = t + (row_number() - 1) * 1e-3) %>%   # break exact same-day ties, negligible vs. day-scale resolution
  ungroup()

N_obs_elife <- nrow(vl_long)

contacts_elife <- df_elife %>%
  left_join(cases_elife %>% select(case_id, case_idx, symp_date), by = "case_id") %>%
  mutate(exp_date = as.Date(date_first_exposure),
         t_exposure = as.numeric(exp_date - symp_date)) %>%
  filter(!is.na(infected_contact), !is.na(t_exposure))

cat(sprintf(
  "eLife contacts: %d total, %d retained after dropping missing exposure/symptom dates (%d dropped)\n",
  nrow(df_elife), nrow(contacts_elife), nrow(df_elife) - nrow(contacts_elife)
))

mydata$N_ind_elife <- N_ind_elife
mydata$N_obs_elife <- N_obs_elife
mydata$y_obs_elife <- log(10^vl_long$vl)   # natural log of linear titer, matches y_obs convention
mydata$ind_elife <- vl_long$case_idx
mydata$t_elife <- vl_long$t
mydata$male_elife <- cases_elife$is.male
mydata$over65_elife <- cases_elife$is.over65

mydata$inf_contact <- contacts_elife$infected_contact
mydata$N_contact <- nrow(contacts_elife)
mydata$contact_case_idx <- contacts_elife$case_idx
mydata$t_exposure_elife <- contacts_elife$t_exposure

mydata$trans_lag_days <- 0    # concurrent-Vt variant: always 0, no CLI override
mydata$d_fixed_elife <- 0.0

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------
out_tag <- "expdecay_concurrent_lag0"
cat(sprintf("Fitting concurrent-Vt expdecay model with trans_lag_days = %s\n", mydata$trans_lag_days))

fit <- stan(
  file = "myclaphammodel_v17_elife_expdecay_concurrent.stan",
  data = mydata,
  seed = 1234565,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  warmup = 4000,
  iter = 5000,
  refresh = 200
)

save(fit, file = sprintf("elife_kinetics_fit_%s.RData", out_tag))
posterior_samples <- rstan::extract(fit)
save(posterior_samples, file = sprintf("elife_kinetics_posterior_%s.RData", out_tag))

summ <- as.data.frame(summary(fit)[[1]])
summ$variable <- row.names(summ)
save(summ, file = sprintf("elife_kinetics_summary_%s.RData", out_tag))

# ---------------------------------------------------------------------------
# Identifiability diagnostics, adapted for the exp-decay parameterization
# (one individual-level parameter, A_elife[i], instead of three).
# ---------------------------------------------------------------------------

## 1. Divergences and Rhat, hyperparameters vs individual-level
diag_elife <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
n_divergent <- sum(sapply(diag_elife, function(x) sum(x[, "divergent__"])))
n_draws_total <- sum(sapply(diag_elife, nrow))
cat(sprintf("\neLife submodel (exp-decay): %d divergent transitions out of %d post-warmup draws\n",
            n_divergent, n_draws_total))

summ_mat <- summary(fit)$summary
is_indiv_elife <- grepl("^A_elife\\[", rownames(summ_mat))
is_hyper_elife <- rownames(summ_mat) %in% c("A_mu_elife", "A_sigma_elife", "k_elife", "sigma_meas_elife")

cat("\nRhat, eLife hyperparameters (population-level, should be clean):\n")
print(summary(summ_mat[is_hyper_elife, "Rhat"]))
cat("\nRhat, eLife individual-level intercepts A_elife[i]:\n")
print(summary(summ_mat[is_indiv_elife, "Rhat"]))
cat(sprintf("Individual-level eLife params with Rhat > 1.01: %d of %d\n",
            sum(summ_mat[is_indiv_elife, "Rhat"] > 1.01, na.rm = TRUE), sum(is_indiv_elife)))

## 2. Prior-vs-posterior contraction, split by how many VL points each case had
n_per_case <- vl_long %>% count(case_idx, name = "n_pts")
all_cases <- tibble(case_idx = 1:N_ind_elife) %>%
  left_join(n_per_case, by = "case_idx") %>%
  mutate(n_pts = replace_na(n_pts, 0))
cat("\neLife cases by number of VL points:\n")
print(table(all_cases$n_pts))

zero_pt_ids <- all_cases %>% filter(n_pts == 0) %>% pull(case_idx) %>% head(5)
one_pt_ids  <- all_cases %>% filter(n_pts == 1) %>% pull(case_idx) %>% head(5)
check_ids_elife <- c(zero_pt_ids, one_pt_ids)

if (length(check_ids_elife) > 0) {
  post_raw <- rstan::extract(fit, pars = "A_raw_elife")$A_raw_elife
  contraction_elife <- map_dfr(check_ids_elife, function(i) {
    tibble(
      case_idx = i,
      n_pts = all_cases$n_pts[all_cases$case_idx == i],
      post_sd_A_raw = sd(post_raw[, i])
    )
  })
  cat("\nPosterior SD of non-centered eLife individual deviations (prior SD = 1;\n",
      "values near 1 -- expected for the 0-point cases especially, since they\n",
      "have no likelihood contribution at all -- mean that individual's intercept\n",
      "is purely the population prior, not case-specific information):\n", sep = "")
  print(contraction_elife)
}

## 3. Pairs plot for one zero-point and one one-point case: A_elife[i] vs the
## shared k_elife, to check for the same kind of intercept/rate tradeoff that
## motivated dropping the up-down curve's f_elife/ci_elife in the first place.
pairs_ids <- c(head(zero_pt_ids, 1), head(one_pt_ids, 1))
for (i in pairs_ids) {
  pars_i <- c(sprintf("A_elife[%d]", i), "k_elife")
  p <- mcmc_pairs(as.array(fit), pars = pars_i, diag_fun = "dens", off_diag_fun = "scatter")
  ggsave(sprintf("elife_kinetics_pairs_case%d_%s.png", i, out_tag), p, width = 8, height = 8)
}
if (length(pairs_ids) > 0) {
  cat("Saved pairs plots for cases:", paste(pairs_ids, collapse = ", "), "\n")
}

## 4. The consequential check: does vt_exposure uncertainty actually shrink
## with more donor-case VL points?
post_vt <- rstan::extract(fit, pars = "vt_exposure")$vt_exposure
vt_sd <- apply(post_vt, 2, sd)
contact_diag <- tibble(
  contact_row = seq_along(vt_sd),
  case_idx = mydata$contact_case_idx,
  vt_exposure_post_sd = vt_sd
) %>%
  left_join(all_cases, by = "case_idx")

cat("\nPosterior SD of vt_exposure (the actual transmission-model covariate), by donor case's # VL points:\n")
contact_diag %>%
  group_by(n_pts) %>%
  summarise(n_contacts = n(), median_post_sd = median(vt_exposure_post_sd), .groups = "drop") %>%
  print()

## 5. The quantity of actual scientific interest: does trans1l (viral-load
## effect on transmission probability) change relative to the up-down fit?
## Compare directly against elife_kinetics_summary_lag<N>.RData if present.
cat("\ntrans0 / trans1l posterior summary (exp-decay eLife submodel):\n")
print(summ_mat[c("trans0", "trans1l"), c("mean", "sd", "2.5%", "97.5%", "Rhat")])

updown_summary_file <- "elife_kinetics_summary_concurrent_lag0.RData"
if (file.exists(updown_summary_file)) {
  e <- new.env()
  load(updown_summary_file, envir = e)
  summ_updown <- e$summ
  cat("\ntrans0 / trans1l posterior summary (up-down eLife submodel, concurrent-Vt, for comparison):\n")
  print(summ_updown[summ_updown$variable %in% c("trans0", "trans1l"),
                     c("variable", "mean", "sd", "2.5%", "97.5%", "Rhat")])
} else {
  cat(sprintf("\n[%s not found -- run the up-down variant with the same trans_lag_days to compare]\n",
              updown_summary_file))
}
