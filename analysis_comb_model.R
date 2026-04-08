library(tidyverse)

dfi <- read_csv2("./dataset_PNAS_SARSCOV2_655_anon_clean.csv")

sum(is.na(dfi$sex))

dfi %>%
  filter(!is.na(sex)) -> dfi_filt

length(unique(dfi_filt$ID))

dfi_filt %>%
  filter(type==1) %>%
  group_by(ID) %>%
  summarise(count = n()) -> dcount

length(dcount$ID)
length(dcount[dcount$count>3,]$ID)

dfi_filt %>%
  filter(type==1) %>%
  group_by(ID) %>%
  arrange(ID, time_monolix) %>%
  mutate(y_lag = lag(y)) %>%
  mutate(diff_y = y - y_lag) %>%
  mutate(miny = min(diff_y, na.rm = T)) %>%
  mutate(maxy = max(diff_y, na.rm = T)) %>%
  mutate(flag = if_else(miny==0 & maxy==0, 1, 0)) %>%
  filter(flag==1) -> dobsequal

sum(dobsequal$flag==1)

sum(dcount$count>3) 
sum(dcount$count<=3) 

dcount %>% filter(count>3) -> xx

length(dcount$ID)
length(unique(xx$ID))

length(unique(dfi_filt$ID))
dfi_filt %>% 
  filter(ID %in% xx$ID) -> dfi_filt_a 
length(unique(dfi_filt_a$ID))

dfi_filt %>% 
  filter(!(ID %in% xx$ID)) -> x99
length(unique(x99$ID))

dfi_filt_a %>%
  filter(ID %in% dobsequal$ID) -> dequal

IDs_equal = unique(dequal$ID)
length(unique(dequal$ID))

length(unique(dfi_filt_a$ID))

dfi_filt_a %>%
  filter(!(ID %in% IDs_equal)) -> dfi_filt_b

dfi_filt_orig  <- dfi_filt

dfi_filt <- dfi_filt_b

dfi_filt %>%
  filter(type==2) %>%
  group_by(ID) %>%
  mutate(died = max(y)) %>%
  filter(time_monolix !=0) %>%
  dplyr::select(ID, died, time_monolix) %>%
  distinct(ID, died, time_monolix) %>%
  mutate(tevent = as.integer(time_monolix)) %>%
  select(!time_monolix) -> df_fatal

table(df_fatal$died)

dfi_filt2 <- left_join(dfi_filt, df_fatal, by = c("ID"="ID"))

length(unique(dfi_filt2$ID))

dfi_filt2 %>%
  filter(type ==1) %>%
  ungroup() %>%
  mutate(vir = 10^y) %>%
  mutate(ID_fct = as.factor(ID)) %>%
  mutate(i_pt = as.numeric(ID_fct)) %>%
  mutate(sero = 1) %>%
#  mutate(is.over65 = if_else(Age>15, 1, 0)) %>%
  mutate(is.over65 = age_cat_reg) %>%
  mutate(is.male = 1-sex) %>%
#  mutate(severe = if_else(str_detect(SevMax3, "Severe"), 1, 0)) %>%
  mutate(severe = died) %>%
#  mutate(severe = if_else(SevMax>=6, 1, 0)) %>%
  group_by(i_pt) %>%
  mutate(Day = time_monolix) %>%
  arrange(Day) %>%
  mutate(Day_lag = lag(Day)) %>%
  mutate(diff_day = if_else(is.na(Day_lag), 1, Day - Day_lag)) %>%
  mutate(prod_vir_dt = vir*diff_day) %>%
  mutate(maxi = max(vir)) %>%
#  filter(maxi >1) %>%
#  filter(StudyNum== 14 | StudyNum== 16 | StudyNum == 10) %>%  # studies with mortality not remdesivir
  mutate(maxv = sum(prod_vir_dt, na.rm=TRUE)) %>%
  mutate(peak_vir = max(vir, na.rm=TRUE)) %>%
  mutate(max_Day = max(Day)) %>%
  mutate(min_Day = min(Day)) %>%
  mutate(days_obs = max_Day - min_Day+1) %>%
  mutate(rate_vir = maxv/days_obs) -> dfc3i_filt

dfc3i_filt %>%
  filter(type==1) -> dfc3_filt

dfc3_filt %>%
  dplyr::select(i_pt, is.male, is.over65, severe, sero, rate_vir, 
                tevent,
                peak_vir, maxv) %>%
  distinct(i_pt, is.male, is.over65, severe, sero, rate_vir, 
           tevent,
           peak_vir, maxv) -> dfu

#Data for popupalion submodel 
# Ile de France

YDF = read_csv("FR_IDF.csv") 

YDFcases <- YDF[70:132, "new_confirmed"]
YDFdeaths <- YDF[70:132, "new_deceased"]

dfYDF <- data.frame(YDFcases = unlist(YDFcases), 
                    YDFdeaths= unlist(YDFdeaths),
                    tc_obs = 1:length(unlist(YDFcases)))

dfYDF %>%
  mutate(cumcases = cumsum(YDFcases)) %>%
  mutate(cumdeaths = cumsum(YDFdeaths)) -> dfYDFcum
  
orig <- 70:132
iorig <- orig -69
iorig
iorig/7

isem = (1:9)*7

ncasescum = dfYDFcum$cumcases[isem]
ndeathscum = dfYDFcum$cumdeaths[isem]

ddiv = 7
dfYDF %>%
  mutate(sem =  1 + (tc_obs-1) %/% ddiv) %>%
  add_count(sem, wt = YDFcases, name = "cases7") %>%
  add_count(sem, wt = YDFdeaths, name = "deaths7") %>%
  dplyr::select(sem, cases7, deaths7) %>%
  distinct(sem, cases7, deaths7) %>%
  mutate(tc_obs = sem*ddiv) -> xxYDF

dfYDF %>%
  filter(tc_obs %% ddiv == 1) -> xx2

xx2$YDFcases

plot(xx2$YDFcases)

plot(unlist(YDF[70:130, "new_confirmed"]))

dfYDF %>%
  replace_na(list(YDFdeaths=0)) -> dfYDF

# if using all 1 by 1
ncases <- unlist(dfYDF$YDFcases)
nd <- unlist(dfYDF$YDFdeaths)
tc_obs <- unlist(dfYDF$tc_obs)

# if aggregating by 7
ncases <- unlist(xxYDF$cases7)
nd <- unlist(xxYDF$deaths7)
tc_obs <- unlist(xxYDF$tc_obs)
tc_obs

npop = unlist(YDF[1, "population"])

mods <- lm(log(cases7) ~ tc_obs, data=xxYDF[1:3,])
coef(mods)[1]
Y0rough_estim=exp(coef(mods)[1])
coef(mods)[2]

gamma_fixed = 1/6.5
CFRrough = sum(nd)/sum(ncases)
mu = CFRrough*gamma_fixed/(1-CFRrough) 

betarough = unname(coef(mods)[2] + mu + gamma_fixed)
betarough

n_t = max(tc_obs)
ts = seq(1, max(tc_obs), by=1) 
t0 = 0
n_tcobs = length(tc_obs)


mydata <- list(y_obs=  log(dfc3_filt$vir),
               ysev = dfu$severe,
               N_ind = max(dfc3_filt$i_pt),
               N_obs = dim(dfc3_filt)[1],
               male = dfu$is.male,
               over65 = dfu$is.over65,
               tevent = dfu$tevent,
               nserotype =1,
               serotype = dfu$sero,
               ratevir = dfu$rate_vir,
               peak = dfu$peak_vir,
               area = dfu$maxv,
               eps0 = -8.91,  # obtained for dengue chaptgpt study denv1 ref. 22 of ben-shachar and koelle, nature communications
               eps1 = 1.368,
               myinf = 30,
               sigma_meas = .5,
               d_fixed = 0.0,
               t = dfc3_filt$Day,
               ind = dfc3_filt$i_pt,
               # parameters for population model
               ts = ts,
               t0 = t0,
              # ncases = log(ncasescum+1.0),
              # ncases = log(ncases+1.0),
                ncases = ncases,
               ndeaths = nd,
 #ndeaths = log(nd+1.0),
 CFRrough = sum(nd)/sum(ncases),
 betarough = betarough,
# ndeaths = log(ndeathscum+1.0),
n_t = n_t,
               tc_obs = tc_obs,
# popnorm = 10000,
popnorm = npop,
npop = npop,
Y0 = Y0rough_estim,
#Y0 = 20,
n_tcobs = n_tcobs
)

dfu %>%
  ggplot(aes(x=log10(maxv), y= severe)) + geom_point() + facet_grid(is.over65 ~ .)

l1 <- glm(severe ~ log10(maxv), family = "binomial", data= dfu)
summary(l1)
confint(l1)
l1 <- glm(severe ~ log10(maxv) + is.male, family = "binomial", data= dfu)
summary(l1)
confint(l1)
l1 <- glm(severe ~ log10(maxv) + is.male + is.over65, family = "binomial", data= dfu)
summary(l1)
confint(l1)
l1 <- glm(severe ~ log10(maxv)*is.male , family = "binomial", data= dfu)
summary(l1)
confint(l1)
l1 <- glm(severe ~ log10(peak_vir ), family = "binomial", data= dfu)
summary(l1)
confint(l1)

#mydata$unobs_idx <- setdiff(1:mydata$N + 1, mydata$obs_idx)

source("./analise_trans_submodel.R")

mydata$inf_contact = df_inf$infected_contact
mydata$day_VL = log(10^df_inf$case_day0_VL)
mydata$N_contact = length(mydata$inf_contact)

library(rstan)
library(ggplot2)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")


myargs <- commandArgs(trailingOnly = TRUE)
print(myargs[1])

modelfile <- 'myclaphammodel_v16_noEIP_covar_surv_ln.stan'
print(modelfile)

nchain <- 4

LOADFIRST = FALSE

mydata$nresp = max(mydata$tevent)

if (!LOADFIRST) {

  adapt = .99
  max_tree = 15
fit <- stan(
  #file = 'myclaphammodel_v6c_EIP.stan',
  file = modelfile,
  data=mydata,
             seed=1234565,
            chains = nchain,
              control = list(
                    adapt_delta = adapt,      # up from default 0.8
                    max_treedepth = max_tree       # up from default 10
                        ),
             warmup = 4000,
              refresh=200,
             # iter = 5000,
             iter = 5000 )

gc()
save(fit, file = "bmc_fit_leakage_ln.RData")
posterior_samples <- rstan::extract(fit)
gc()

#fit

summ1 <- summary(fit)
summ <- summ1[[1]]

summdf <- as.data.frame(summ)
summdf$variable = row.names(summdf)

save(summdf, file = "bmc_leakage_ln.RData")
save(posterior_samples, file = "bmc_leakage_posterior_ln.RData")
}
#load("bmc_fit_leakage.RData")
#load("bmc_leakage_ln.RData")
#load("bmc_fit_leakage_ln.RData")
load("bmc_leakage_posterior_ln.RData")

posterior <- posterior_samples

scaled_sumrv <- sweep(posterior$sumrv, 1, posterior$constant_R0, `*`)

n_iter <- nrow(scaled_sumrv)
n_ind  <- ncol(scaled_sumrv)

df_all <- tibble(
  id    = rep(seq_len(n_ind), each = n_iter),
  iter  = rep(seq_len(n_iter), times = n_ind),
  psev  = as.vector(posterior$psev),       # [iter x 81] -> vector, column-major
  sumrv = as.vector(scaled_sumrv)
)

library(patchwork)

df_rect <- tibble(
  id         = seq_len(n_ind),
  psev_median = apply(posterior$psev, 2, median),
  psev_lo95  = apply(posterior$psev,  2, quantile, 0.025),
  psev_hi95  = apply(posterior$psev,  2, quantile, 0.975),
  psev_lo50  = apply(posterior$psev,  2, quantile, 0.25),
  psev_hi50  = apply(posterior$psev,  2, quantile, 0.75),
  sumrv_median = apply(scaled_sumrv,  2, median),
  sumrv_lo95 = apply(scaled_sumrv,    2, quantile, 0.025),
  sumrv_hi95 = apply(scaled_sumrv,    2, quantile, 0.975),
  sumrv_lo50 = apply(scaled_sumrv,    2, quantile, 0.25),
  sumrv_hi50 = apply(scaled_sumrv,    2, quantile, 0.75)
)

df_meta <- tibble(
  id     = seq_len(n_ind),
  male   = mydata$male,
  over65 = mydata$over65
) %>%
  mutate(
    sex_label = if_else(male   == 1, "Male",  "Female"),
    age_label = if_else(over65 == 1, "≥65 yrs", "<65 yrs"),
    group     = paste(sex_label, age_label, sep = ", ")
  )

df_rect <- df_rect %>% left_join(df_meta, by = "id")

df_all  <- df_all  %>%
  left_join(df_meta, by = c("id" = "id"))

group_colours <- c(
  "Male, ≥65 yrs"   = "#2166ac",
  "Male, <65 yrs"   = "#92c5de",
  "Female, ≥65 yrs" = "#b2182b",
  "Female, <65 yrs" = "#f4a582"
)

group_shapes <- c(
  "Male, ≥65 yrs"   = 16,
  "Male, <65 yrs"   = 17,
  "Female, ≥65 yrs" = 15,
  "Female, <65 yrs" = 18
)

facet_theme <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    plot.margin      = margin(8, 8, 8, 8)
  )

##

group_shapes <- c(
  "Male, ≥65 yrs"   = 16,
  "Male, <65 yrs"   = 17,
  "Female, ≥65 yrs" = 15,
  "Female, <65 yrs" = 18
)

p_scatter <- ggplot() +
  geom_point(
    data  = df_all,
    aes(x = psev, y = sumrv, colour = group),
    size  = 0.2, alpha = 0.02
  ) +
  geom_density_2d(
    data      = df_all,
    aes(x = psev, y = sumrv, colour = group),
    linewidth = 0.4,
    alpha     = 0.9
  ) +
  scale_colour_manual(
    values = group_colours,
    name   = NULL,
    guide  = guide_legend(
      override.aes   = list(size = 3, alpha = 1,
                         #   shape = group_shapes,
                            linetype = 0),
      title.position = "top",
      nrow           = 2
    )
  ) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 5),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylim(c(0, 5)) +
  labs(
    x     = expression(italic(FP)),
    y     = expression(pR[0]),
    title = "Posterior samples"
  )

p_box_facet <- ggplot() +
  geom_rect(
    data = df_rect,
    aes(xmin = psev_lo95, xmax = psev_hi95,
        ymin = sumrv_lo95, ymax = sumrv_hi95,
        colour = group),
    fill = NA, linewidth = 0.2
  ) +
  geom_rect(
    data = df_rect,
    aes(xmin = psev_lo50, xmax = psev_hi50,
        ymin = sumrv_lo50, ymax = sumrv_hi50,
        fill = group, colour = group),
    linewidth = 0.3, alpha = 0.25
  ) +
  geom_point(
    data  = df_rect,
    aes(x = psev_median, y = sumrv_median),
    shape = 3, size = 1.2, stroke = 0.6, color = "black"
  ) +
  scale_colour_manual(values = group_colours, name = NULL) +
  scale_fill_manual(  values = group_colours, name = NULL, guide = "none") +
  facet_grid(age_label ~ sex_label) +
  scale_x_log10(
    breaks = scales::log_breaks(n = 5),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  ylim(c(0, 5)) +
  labs(
    x= "FP",
    y = expression(italic(pR)[0]),
  #  x     = expression(italic(p)[sev]),
  #  y     = expression(beta * Sigma*italic(rv)),
    title = "50% and 95% credible regions"
  )

my_theme_pub <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    plot.margin      = margin(8, 8, 8, 8)
  )

my_theme_pub_x <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "none",
    legend.text      = element_text(size = 9),
    plot.margin      = margin(8, 8, 8, 8)
  )

p_combined1 <- (p_scatter + my_theme_pub) +
  (p_box_facet + my_theme_pub_x) +
  plot_layout(guides = "collect", widths = c(1, 1.4)) +
  plot_annotation(
   # caption = "Left: all posterior draws coloured by group. Right: 50% (filled) and 95% (outline) credible regions faceted by sex and age.",
    theme   = theme(
      plot.caption    = element_text(size = 8, colour = "grey50", hjust = 0),
      legend.position = "bottom"
    )
  )

p_combined1

my_theme_pub <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 11),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.title     = element_blank(),
    plot.margin      = margin(8, 8, 8, 8)
  )

ggsave("psev_sumrv_combined.pdf",  plot = p_combined1, width = 11, height = 6, device = cairo_pdf)
ggsave("psev_sumrv_combined.tiff", plot = p_combined1, width = 11, height = 6, dpi = 600, compression = "lzw")


library(gt)

df_summary_table <- df_all %>%
  group_by(group, sex_label, age_label) %>%
  summarise(
    # scaled sumrv
    sumrv_median = median(sumrv,          na.rm = TRUE),
    sumrv_q25    = quantile(sumrv, 0.25,  na.rm = TRUE),
    sumrv_q75    = quantile(sumrv, 0.75,  na.rm = TRUE),
    sumrv_lo95   = quantile(sumrv, 0.025, na.rm = TRUE),
    sumrv_hi95   = quantile(sumrv, 0.975, na.rm = TRUE),
    # psev
    psev_median  = median(psev,           na.rm = TRUE),
    psev_q25     = quantile(psev, 0.25,   na.rm = TRUE),
    psev_q75     = quantile(psev, 0.75,   na.rm = TRUE),
    psev_lo95    = quantile(psev, 0.025,  na.rm = TRUE),
    psev_hi95    = quantile(psev, 0.975,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  mutate(
    sumrv_iqr  = glue::glue("[{sumrv_q25}, {sumrv_q75}]"),
    sumrv_ci95 = glue::glue("[{sumrv_lo95}, {sumrv_hi95}]"),
    psev_iqr   = glue::glue("[{psev_q25}, {psev_q75}]"),
    psev_ci95  = glue::glue("[{psev_lo95}, {psev_hi95}]")
  ) %>%
  select(group, sex_label, age_label,
         sumrv_median, sumrv_iqr, sumrv_ci95,
         psev_median,  psev_iqr,  psev_ci95)

tbl_summary <- df_summary_table %>%
  select(-sex_label, -age_label) %>%
  gt(rowname_col = "group") %>%
  tab_spanner(
    label   = md("**β · Σrv** (scaled sumrv)"),
    columns = c(sumrv_median, sumrv_iqr, sumrv_ci95)
  ) %>%
  tab_spanner(
    label   = md("***p*_sev**"),
    columns = c(psev_median, psev_iqr, psev_ci95)
  ) %>%
  cols_label(
    sumrv_median = "Median",
    sumrv_iqr    = "IQR [25%, 75%]",
    sumrv_ci95   = "95% interval",
    psev_median  = "Median",
    psev_iqr     = "IQR [25%, 75%]",
    psev_ci95    = "95% interval"
  ) %>%
  tab_header(
    title    = md("**Summary statistics by sex and age group**"),
    subtitle = md("Posterior samples of β·Σrv and *p*_sev")
  ) %>%
  # colour row stubs to match plot colours
  tab_style(
    style     = list(cell_fill(color = "#2166ac"),
                     cell_text(color = "white", weight = "bold")),
    locations = cells_stub(rows = group == "Male, ≥15 yrs")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#92c5de"),
                     cell_text(color = "white", weight = "bold")),
    locations = cells_stub(rows = group == "Male, <15 yrs")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#b2182b"),
                     cell_text(color = "white", weight = "bold")),
    locations = cells_stub(rows = group == "Female, ≥15 yrs")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#f4a582"),
                     cell_text(color = "#333333", weight = "bold")),  # dark text for light background
    locations = cells_stub(rows = group == "Female, <15 yrs")
  ) %>%
  tab_style(
    style     = cell_fill(color = "#f0f4ff"),
    locations = cells_body(columns = c(sumrv_median, sumrv_iqr, sumrv_ci95))
  ) %>%
  tab_style(
    style     = cell_fill(color = "#fff4f0"),
    locations = cells_body(columns = c(psev_median, psev_iqr, psev_ci95))
  ) %>%
  tab_footnote(
    footnote  = "Statistics computed across all posterior samples within each group.",
    locations = cells_title("subtitle")
  ) %>%
  cols_align(align = "center",
             columns = c(sumrv_median, sumrv_iqr, sumrv_ci95,
                         psev_median,  psev_iqr,  psev_ci95)) %>%
  cols_align(align = "left", columns = everything()) %>%
  tab_options(
    table.font.names              = "serif",
    table.font.size               = 12,
    heading.align                 = "left",
    column_labels.border.top.width    = px(2),
    column_labels.border.bottom.width = px(2),
    table.border.top.width            = px(2),
    table.border.bottom.width         = px(2),
    data_row.padding                  = px(6),
    stub.border.width                 = px(0)
  )

gtsave(tbl_summary, "sumrv_psev_summary_table.html")
gtsave(tbl_summary, "sumrv_psev_summary_table.pdf")
gtsave(tbl_summary, "sumrv_psev_summary_table.rtf")

print(tbl_summary)


n_samples <- nrow(posterior$psev)
n_ind     <- ncol(posterior$psev)

df_psev_maxi <- tibble(
  individual = rep(seq_len(n_ind), each = n_samples),
  sample     = rep(seq_len(n_samples), times = n_ind),
  psev       = as.vector(posterior$psev),
  maxi       = as.vector(posterior$maxi)
)

my_theme <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 10),
    plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    plot.margin      = margin(8, 8, 8, 8)
  )

group_colours <- c(
  "Male, ≥65 yrs"   = "#2166ac",
  "Male, <65 yrs"   = "#92c5de",
  "Female, ≥65 yrs" = "#b2182b",
  "Female, <65 yrs" = "#f4a582"
)

library(deSolve)

n_samples <- nrow(posterior$beta_ind)   # number of posterior draws
n_ind     <- ncol(posterior$beta_ind)   # number of individuals

# Individual-level parameters [n_samples x n_ind]
beta_ind  <- posterior$beta_ind
gamma_ind <- posterior$gamma_ind
omega_ind <- posterior$omega_ind
pop_ind   <- posterior$pop_ind

# Scalar across individuals [n_samples]
Y0f <- posterior$Y0f

sir_ode <- function(t, state, params) {
  S  <- state["S"]
  I  <- state["I"]
  R  <- state["R"]
  Cd <- state["Cd"]   # cumulative deaths
  Cc <- state["Cc"]   # cumulative cases
  
  beta  <- params["beta"]
  gamma <- params["gamma"]
  omega <- params["omega"]
  pop   <- params["pop"]
  
  N <- S + I + R
  dS  <- -beta * S * I / N
  dI  <-  beta * S * I / N - gamma * I - omega * I
  dR  <-  gamma * I
  dCd <-  omega * I
  dCc <-  beta  * S * I / N
  
  list(c(dS, dI, dR, dCd, dCc))
}

my_tf = 63

##

threshold <- quantile(df_all$sumrv, 0.90)

df_rect <- df_rect %>%
 # rename(individual = id) %>%
  left_join(
    df_all %>%
      group_by(id) %>%
      summarise(
        p_exceed  = mean(sumrv > threshold),   # posterior prob of exceeding
        sumrv_med = median(sumrv),
        .groups   = "drop"
      ),
    by = "id"
  )

df_rect <- df_rect %>%
  left_join(
    df_all %>%
      group_by(id) %>%
      summarise(
        sumrv_med = median(sumrv),
        # density at the median = concentration of posterior mass
        dens_at_med = {
          d <- density(sumrv, n = 512)
          approx(d$x, d$y, xout = median(sumrv))$y
        },
        .groups = "drop"
      ),
    by = "id"
  )

df_rect <- df_rect %>%
  mutate(
    sumrv_score = sumrv_median * p_exceed   # high median AND high certainty
  )


### observations population level


y_rep  <- posterior$y_rep    # [n_samples x length(tc_obs)]
yd_rep <- posterior$yd_rep   # [n_samples x length(tc_obs)]

df_obs <- tibble(
  time   = tc_obs,
  cases  = mydata$ncases,
  deaths = mydata$ndeaths
)

my_theme <- theme_classic(base_size = 12, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 10),
    axis.title       = element_text(size = 12),
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    plot.margin      = margin(10, 10, 10, 10)
  )

df_fit <- tibble(
  time          = tc_obs,
  cases_median  = apply(y_rep,  2, median),
  cases_lo50    = apply(y_rep,  2, quantile, 0.25),
  cases_hi50    = apply(y_rep,  2, quantile, 0.75),
  cases_lo95    = apply(y_rep,  2, quantile, 0.025),
  cases_hi95    = apply(y_rep,  2, quantile, 0.975),
  deaths_median = apply(yd_rep, 2, median),
  deaths_lo50   = apply(yd_rep, 2, quantile, 0.25),
  deaths_hi50   = apply(yd_rep, 2, quantile, 0.75),
  deaths_lo95   = apply(yd_rep, 2, quantile, 0.025),
  deaths_hi95   = apply(yd_rep, 2, quantile, 0.975)
)

plot_fit <- function(obs_var, median_var,
                     lo50_var, hi50_var,
                     lo95_var, hi95_var,
                     y_label, title_label) {
  ggplot() +
    # 95% CrI — outer, lighter
    geom_ribbon(
      data = df_fit,
      aes(x = time, ymin = .data[[lo95_var]], ymax = .data[[hi95_var]]),
      fill = "#2166ac", alpha = 0.15
    ) +
    # 50% IQR — inner, darker
    geom_ribbon(
      data = df_fit,
      aes(x = time, ymin = .data[[lo50_var]], ymax = .data[[hi50_var]]),
      fill = "#2166ac", alpha = 0.35
    ) +
    # Posterior median line
    geom_line(
      data = df_fit,
      aes(x = time, y = .data[[median_var]]),
      colour = "#2166ac", linewidth = 0.7
    ) +
    # Observed points
    geom_point(
      data = df_obs,
      aes(x = time, y = .data[[obs_var]]),
      shape = 21, size = 1.8,
      fill = "white", colour = "#b2182b", stroke = 0.5
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
    scale_y_continuous(labels = scales::comma) +
    labs(x = "Days", y = y_label, title = title_label) +
    my_theme +
    theme(plot.title = element_text(size=12, hjust = -0.2))
}

p_cases_fit <- plot_fit(
  obs_var     = "cases",
  median_var  = "cases_median",
  lo50_var    = "cases_lo50",
  hi50_var    = "cases_hi50",
  lo95_var    = "cases_lo95",
  hi95_var    = "cases_hi95",
  y_label     = "Number of cases",
  title_label = "A"
)

p_deaths_fit <- plot_fit(
  obs_var     = "deaths",
  median_var  = "deaths_median",
  lo50_var    = "deaths_lo50",
  hi50_var    = "deaths_hi50",
  lo95_var    = "deaths_lo95",
  hi95_var    = "deaths_hi95",
  y_label     = "Number of deaths",
  title_label = "B"
)

p_fit <- p_cases_fit + p_deaths_fit +
  plot_annotation(
  #  caption = "Blue line: posterior median. Dark band: 50% IQR. Light band: 95% credible interval. Red circles: observed data.",
    theme   = theme(plot.caption = element_text(size = 8, colour = "grey50", hjust = 0))
  )

p_fit 

ggsave("observed_vs_fit.pdf",  plot = p_fit, width = 10, height = 4.5, device = cairo_pdf)
ggsave("observed_vs_fit.tiff", plot = p_fit, width = 10, height = 4.5, dpi = 600, compression = "lzw")


# table #

#library(rstan)
#library(tidyverse)
#library(gt)

# Format a vector of samples as "median [lo, hi]"
fmt_cri <- function(x, digits = 2, uselog10=FALSE) {
  if (uselog10) {
    med <- round(log(10)*median(x, na.rm = TRUE), digits)
    lo  <- round(log(10)*quantile(x, 0.025, na.rm = TRUE), digits)
    hi  <- round(log(10)*quantile(x, 0.975, na.rm = TRUE), digits)
  } else {
  med <- round(median(x, na.rm = TRUE), digits)
  lo  <- round(quantile(x, 0.025, na.rm = TRUE), digits)
  hi  <- round(quantile(x, 0.975, na.rm = TRUE), digits)
  }
  glue::glue("{med} [{lo}, {hi}]")
}

# For individual-level: median per individual, then mean of medians and range
fmt_ind <- function(mat, digits = 2) {
  # mat is [n_samples x n_individuals]
  medians <- apply(mat, 2, median, na.rm = TRUE)
  m       <- round(mean(medians), digits)
  lo      <- round(min(medians),  digits)
  hi      <- round(max(medians),  digits)
  glue::glue("{m} [{lo}, {hi}]")
}


df_table <- tribble(
  ~Model,           ~Parameter,                                        ~Summary,
  
  # --- Transmission model ---
  "Transmission",   "trans0",                                          fmt_cri(posterior$trans0),
  "Transmission",   "trans1",                                          fmt_cri(posterior$trans1, uselog10 = TRUE),
  
  # --- Survival model ---
  "Survival",       "gamma0 (intercept)",                              fmt_cri(posterior$gamma0),
  "Survival",       "gammaa (sex: male)",                              fmt_cri(posterior$gammaa),
  "Survival",       "gammab (age level binary)",                       fmt_cri(posterior$gammab),
  "Survival",       "gammat (time)",                 fmt_cri(posterior$betat),
  "Survival",       "gammar (cumulative viral titer)",                 fmt_cri(posterior$gammar, uselog10 = TRUE),
  
  # --- Kinetics model ---
  "Kinetics",       "alphas (mean of medians [range])",                fmt_ind(posterior$alphas),
  "Kinetics",       "betas (mean of medians [range])",                 fmt_ind(posterior$betas),
  "Kinetics",       "beta (baseline log-scale)",                       fmt_cri(posterior$beta),
  "Kinetics",       "psisex (sex effect: male)",                       fmt_cri(posterior$psisex),
  "Kinetics",       "psiage (age effect: over 15)",                    fmt_cri(posterior$psiage),
  
  # --- Population-level model ---
  "Population",     "beta_trans",                                      fmt_cri(posterior$beta_trans),
  "Population",     "gamma_trans",                                     fmt_cri(posterior$gamma_trans),
  "Population",     "omega_trans",                                     fmt_cri(posterior$omega_trans),
  "Population",     "R0_trans",                                        fmt_cri(posterior$R0_trans)
)

tbl <- df_table %>%
  gt(groupname_col = "Model") %>%
  cols_label(
    Parameter = "Parameter",
    Summary   = "Median [95% CrI]"
  ) %>%
  tab_header(
    title    = "Posterior parameter estimates",
    subtitle = "Median and 95% credible intervals"
  ) %>%
  tab_footnote(
    footnote = "Kinetics parameters: mean of posterior medians across individuals. Range in brackets.",
    locations = cells_row_groups("Kinetics")
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style     = cell_text(style = "italic"),
    locations = cells_body(columns = Parameter)
  ) %>%
  tab_options(
    table.font.names        = "serif",
    table.font.size         = 12,
    heading.align           = "left",
    row_group.border.top.width    = px(2),
    row_group.border.bottom.width = px(1),
    table.border.top.width        = px(2),
    table.border.bottom.width     = px(2)
  ) %>%
  cols_align(align = "left",  columns = Parameter) %>%
  cols_align(align = "right", columns = Summary)

# As HTML (interactive, high quality)
gtsave(tbl, "posterior_table.html")

# As PDF via webshot2
# install.packages("webshot2")
gtsave(tbl, "posterior_table.pdf")

# As Word-ready RTF
gtsave(tbl, "posterior_table.rtf")

print(tbl)

### kinetics

library(rstan)
library(tidyverse)
library(glue)

#posterior <- rstan::extract(fit)
posterior <- posterior_samples

vt     <- posterior$vt
n_ind  <- dim(vt)[2]
n_resp <- dim(vt)[3]

vt_median <- apply(vt, c(2, 3), median)
vt_lo50   <- apply(vt, c(2, 3), quantile, 0.25)
vt_hi50   <- apply(vt, c(2, 3), quantile, 0.75)
vt_lo95   <- apply(vt, c(2, 3), quantile, 0.025)
vt_hi95   <- apply(vt, c(2, 3), quantile, 0.975)

# Posterior predictive: add measurement noise analytically
pred_lo95 <- vt_lo95 - 1.96 * mydata$sigma_meas
pred_hi95 <- vt_hi95 + 1.96 * mydata$sigma_meas

t_grid <- seq_len(n_resp)   # 1:30

library(glue)

pivot_matrix <- function(mat, value_name) {
  as_tibble(mat) %>%
    mutate(individual = seq_len(n_ind)) %>%
    pivot_longer(-individual,
                 names_to  = "timepoint",
                 values_to = value_name) %>%
    mutate(timepoint = as.integer(str_remove(timepoint, "V")))
}

df_vt <- pivot_matrix(vt_median, "vt_median") %>%
  left_join(pivot_matrix(vt_lo50,   "vt_lo50"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_hi50,   "vt_hi50"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_lo95,   "vt_lo95"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_hi95,   "vt_hi95"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(pred_lo95, "pred_lo95"), by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(pred_hi95, "pred_hi95"), by = c("individual", "timepoint")) %>%
  mutate(t = t_grid[timepoint])


df_obs <- tibble(
  individual = mydata$ind,
  t          = mydata$t,
  y_obs      = mydata$y_obs,
  male       = mydata$male[mydata$ind],
  over65     = mydata$over65[mydata$ind]
) %>%
  mutate(
    sex_label   = if_else(male   == 1, "M", "F"),
    age_label   = if_else(over65 == 1, "≥65", "<65"),
  #  facet_label = glue("Ind {individual} ({sex_label}, {age_label})")
    facet_label = glue("Ind {individual}")
  )


df_vt <- df_vt %>%
  left_join(
    df_obs %>% distinct(individual, facet_label),
    by = "individual"
  )

# For individuals with no observations, fill facet_label manually
df_vt <- df_vt %>%
  mutate(facet_label = if_else(
    is.na(facet_label),
    glue("Ind {individual}"),
    facet_label
  ))

my_theme_small <- theme_classic(base_size = 7, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.2),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(size = 10),
    strip.text       = element_text(size = 11, face = "bold"),
    strip.background = element_blank(),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.2),
    plot.margin      = margin(2, 2, 2, 2),
    plot.caption     = element_text(size = 7, colour = "grey50", hjust = 0)
  )

p_kinetics <- ggplot() +
  # Posterior predictive band (lightest)
  geom_ribbon(
    data = df_vt,
    aes(x = t, ymin = pred_lo95/log(10), ymax = pred_hi95/log(10)),
    fill = "#2166ac", alpha = 0.10
  ) +
  # 95% latent CrI
  geom_ribbon(
    data = df_vt,
    aes(x = t, ymin = vt_lo95/log(10), ymax = vt_hi95/log(10)),
    fill = "#2166ac", alpha = 0.20
  ) +
  # 50% latent IQR
  geom_ribbon(
    data = df_vt,
    aes(x = t, ymin = vt_lo50/log(10), ymax = vt_hi50/log(10)),
    fill = "#2166ac", alpha = 0.35
  ) +
  # Posterior median line
  geom_line(
    data = df_vt,
    aes(x = t, y = vt_median/log(10)),
    colour = "#2166ac", linewidth = 0.5
  ) +
  # Observed points
  geom_point(
    data = df_obs,
    aes(x = t, y = y_obs/log(10)),
    shape = 21, size = 1.5,
    fill = "white", colour = "#b2182b", stroke = 0.4
  ) +
  facet_wrap(
    ~ facet_label,
    nrow   = 9,
    ncol   = 9,
    scales = "free_y"
  ) +
  coord_cartesian(ylim = c(0, NA))+
  labs(
    x       = "Time",
    y       = "Viral titer"
  #  caption = "Blue line: posterior median. Bands: 50% IQR, 95% CrI, 95% predictive interval. Red circles: observed data."
  ) +
  my_theme_small

p_kinetics

ggsave("kinetics_individuals.pdf",
       plot   = p_kinetics,
       width  = 16,
       height = 16,
       device = cairo_pdf)

ggsave("kinetics_individuals.tiff",
       plot        = p_kinetics,
       width       = 16,
       height      = 16,
       dpi         = 300,
       compression = "lzw")

# table 1

fmt_median_iqr <- function(x, digits = 1) {
  med <- round(median(x, na.rm = TRUE), digits)
  lo  <- round(quantile(x, 0.25, na.rm = TRUE), digits)
  hi  <- round(quantile(x, 0.75, na.rm = TRUE), digits)
  glue::glue("{med} [{lo}, {hi}]")
}

fmt_n_pct <- function(x, total) {
  n   <- sum(x, na.rm = TRUE)
  pct <- round(100 * n / total, 1)
  glue::glue("{n} ({pct}%)")
}

# Stratified wrapper: returns overall | male | female
summarise_cat <- function(x, male_vec) {
  n_tot  <- length(x)
  n_m    <- sum(male_vec == 1)
  n_f    <- sum(male_vec == 0)
  list(
    overall = fmt_n_pct(x,              n_tot),
    male    = fmt_n_pct(x[male_vec==1], n_m),
    female  = fmt_n_pct(x[male_vec==0], n_f)
  )
}

summarise_cont <- function(x, male_vec, digits = 1) {
  list(
    overall = fmt_median_iqr(x,              digits),
    male    = fmt_median_iqr(x[male_vec==1], digits),
    female  = fmt_median_iqr(x[male_vec==0], digits)
  )
}

# Individual-level frame
df_ind <- tibble(
  id     = seq_len(mydata$N_ind),
  male   = mydata$male,
  over65 = mydata$over65
)

# Number of observations per individual
n_obs_per_ind <- as.integer(table(mydata$ind))

# Peak observed viral titer per individual
peak_obs <- tapply(mydata$y_obs/log(10), mydata$ind, max, na.rm = TRUE)

kin_male <- df_ind$male

N_kin   <- nrow(df_ind)
N_m_kin <- sum(kin_male == 1)
N_f_kin <- sum(kin_male == 0)

# ysev: 1 = death, 0 = censored (length N_ind)
# tevent: time of death or censoring (length N_ind)

sev_cat  <- summarise_cat(mydata$ysev,    kin_male)
time_evt <- summarise_cont(mydata$tevent, kin_male)
obs_n    <- summarise_cont(n_obs_per_ind, kin_male, digits = 0)
peak_vt  <- summarise_cont(peak_obs,      kin_male)

# inf_contact: binary — transmission occurred (length = n contacts)
# day_VL:      viral titer at day of contact
# These are NOT indexed by individual — summarise overall + by male

# Attempt to get individual index for transmission observations
if (!is.null(mydata$ind_trans)) {
  male_trans <- mydata$male[mydata$ind_trans]
} else {
  male_trans <- rep(NA, length(mydata$inf_contact))
}

trans_cat  <- if (any(!is.na(male_trans))) {
  summarise_cat(mydata$inf_contact, male_trans)
} else {
  list(
    overall = fmt_n_pct(mydata$inf_contact, length(mydata$inf_contact)),
    male    = "—",
    female  = "—"
  )
}

dayvl_cont <- if (any(!is.na(male_trans))) {
  summarise_cont(mydata$day_VL/log(10), male_trans)
} else {
  list(
    overall = fmt_median_iqr(mydata$day_VL/log(10)),
    male    = "—",
    female  = "—"
  )
}

N_trans <- length(mydata$inf_contact)

# tc_obs:    time points (days)
# Cumulative cases and deaths at each tc_obs
# Summarise as median [IQR] across time points
# (no sex stratification — population level)

pop_cases  <- fmt_median_iqr(mydata$ncases)
pop_deaths <- fmt_median_iqr(mydata$ndeaths)
pop_tobs   <- glue::glue("{min(mydata$tc_obs)}–{max(mydata$tc_obs)}")
pop_nobs   <- length(mydata$tc_obs)

rows <- tribble(
  ~section,       ~variable,                              ~overall,                                  ~male,                             ~female,
  
  # ── Kinetics & Survival ──────────────────────────────────────────────────
  "Kinetics & Survival Cohort",
  glue::glue("N individuals"),            as.character(N_kin),                       as.character(N_m_kin),             as.character(N_f_kin),
  "Kinetics & Survival Cohort",
  "Age ≥65 years, n (%)",                 as.character(summarise_cat(df_ind$over65, kin_male)$overall),
  as.character(summarise_cat(df_ind$over65, kin_male)$male),
  as.character(summarise_cat(df_ind$over65, kin_male)$female),
  "Kinetics & Survival Cohort",
  "Observations per individual, median [IQR]",
  as.character(obs_n$overall),               as.character(obs_n$male),          as.character(obs_n$female),
  "Kinetics & Survival Cohort",
  "Peak observed viral titer, median [IQR]",
  as.character(peak_vt$overall),             as.character(peak_vt$male),        as.character(peak_vt$female),
  "Kinetics & Survival Cohort",
  "Death (severe outcome), n (%)",        as.character(sev_cat$overall),             as.character(sev_cat$male),        as.character(sev_cat$female),
  "Kinetics & Survival Cohort",
  "Time to event, days, median [IQR]",    as.character(time_evt$overall),            as.character(time_evt$male),       as.character(time_evt$female),
  
  # ── Transmission ─────────────────────────────────────────────────────────
  "Transmission Cohort",
  glue::glue("N contact events"),         as.character(N_trans),                     "—",                               "—",
  "Transmission Cohort",
  "Transmission occurred, n (%)",         as.character(trans_cat$overall),           as.character(trans_cat$male),      as.character(trans_cat$female),
  "Transmission Cohort",
  "Viral titer at contact, median [IQR]", as.character(dayvl_cont$overall),          as.character(dayvl_cont$male),     as.character(dayvl_cont$female),
  
  # ── Population ───────────────────────────────────────────────────────────
  "Population Level (Île-de-France)",
  "Observation time points (days)",       as.character(pop_tobs),                    "—",                               "—",
  "Population Level (Île-de-France)",
  "N time points",                        as.character(pop_nobs),                    "—",                               "—",
  "Population Level (Île-de-France)",
  "Cumulative cases per interval, median [IQR]",
  as.character(pop_cases),                   "—",                               "—",
  "Population Level (Île-de-France)",
  "Cumulative deaths per interval, median [IQR]",
  as.character(pop_deaths),                  "—",                               "—"
)

tbl <- rows %>%
  gt(groupname_col = "section") %>%
  cols_label(
    variable = "Variable",
    overall  = glue::glue("Overall (N = {N_kin})"),
    male     = glue::glue("Male (N = {N_m_kin})"),
    female   = glue::glue("Female (N = {N_f_kin})")
  ) %>%
  tab_header(
    title    = "Table 1. Characteristics of study cohorts",
    subtitle = "Continuous variables: median [IQR]. Categorical: n (%)."
  ) %>%
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style     = cell_fill(color = "#f7f7f7"),
    locations = cells_row_groups()
  ) %>%
  tab_style(
    style     = cell_text(style = "italic"),
    locations = cells_body(
      columns = variable,
      rows    = str_detect(variable, "N ")
    )
  ) %>%
  tab_footnote(
    footnote  = "Kinetics and survival data come from the same cohort of individuals.",
    locations = cells_row_groups("Kinetics & Survival Cohort")
  ) %>%
  tab_footnote(
    footnote  = "Sex stratification not available for population-level data.",
    locations = cells_row_groups("Population Level (Île-de-France)")
  ) %>%
  cols_align(align = "left",   columns = variable) %>%
  cols_align(align = "center", columns = c(overall, male, female)) %>%
  tab_options(
    table.font.names              = "serif",
    table.font.size               = 12,
    heading.align                 = "left",
    row_group.border.top.width    = px(2),
    row_group.border.bottom.width = px(1),
    table.border.top.width        = px(2),
    table.border.bottom.width     = px(2),
    data_row.padding              = px(5)
  )

gtsave(tbl, "table1.html")
gtsave(tbl, "table1.pdf")
gtsave(tbl, "table1.rtf")   # Word-ready

print(tbl)

##

posterior <- posterior_samples

vt     <- posterior$vt
n_ind  <- dim(vt)[2]
n_resp <- dim(vt)[3]

vt_median <- apply(vt, c(2, 3), median)
vt_lo50   <- apply(vt, c(2, 3), quantile, 0.25)
vt_hi50   <- apply(vt, c(2, 3), quantile, 0.75)
vt_lo95   <- apply(vt, c(2, 3), quantile, 0.025)
vt_hi95   <- apply(vt, c(2, 3), quantile, 0.975)

# Posterior predictive: add measurement noise analytically
pred_lo95 <- vt_lo95 - 1.96 * mydata$sigma_meas
pred_hi95 <- vt_hi95 + 1.96 * mydata$sigma_meas

t_grid <- seq_len(n_resp)   # 1:30

pivot_matrix <- function(mat, value_name) {
  as_tibble(mat) %>%
    mutate(individual = seq_len(n_ind)) %>%
    pivot_longer(-individual,
                 names_to  = "timepoint",
                 values_to = value_name) %>%
    mutate(timepoint = as.integer(str_remove(timepoint, "V")))
}

df_vt <- pivot_matrix(vt_median, "vt_median") %>%
  left_join(pivot_matrix(vt_lo50,   "vt_lo50"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_hi50,   "vt_hi50"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_lo95,   "vt_lo95"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(vt_hi95,   "vt_hi95"),   by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(pred_lo95, "pred_lo95"), by = c("individual", "timepoint")) %>%
  left_join(pivot_matrix(pred_hi95, "pred_hi95"), by = c("individual", "timepoint")) %>%
  mutate(t = t_grid[timepoint])

ind_sumrv_low  <- df_rect %>% slice_min(sumrv_lo50,  n = 1) %>% pull(id)
ind_sumrv_high <- df_rect %>% slice_max(sumrv_hi50,  n = 1) %>% pull(id)
ind_psev_low   <- df_rect %>% slice_min(psev_lo50,   n = 1) %>% pull(id)
ind_psev_high  <- df_rect %>% slice_max(psev_hi50,   n = 1) %>% pull(id)

cat("sumrv IQR low  — individual:", ind_sumrv_low,
    " (Q25 =", round(df_rect$sumrv_lo50[ind_sumrv_low], 3), ")\n")
cat("sumrv IQR high — individual:", ind_sumrv_high,
    " (Q75 =", round(df_rect$sumrv_hi50[ind_sumrv_high], 3), ")\n")
cat("psev  IQR low  — individual:", ind_psev_low,
    " (Q25 =", round(df_rect$psev_lo50[ind_psev_low], 4), ")\n")
cat("psev  IQR high — individual:", ind_psev_high,
    " (Q75 =", round(df_rect$psev_hi50[ind_psev_high], 4), ")\n")

# ============================================================
# Find the sample closest to the IQR boundary for each individual
# ============================================================

# sumrv_low: sample closest to Q25 of individual ind_sumrv_low
target_sumrv_low  <- df_rect$sumrv_lo50[df_rect$id == ind_sumrv_low]
target_sumrv_high <- df_rect$sumrv_hi50[df_rect$id == ind_sumrv_high]
target_psev_low   <- df_rect$psev_lo50[ df_rect$id == ind_psev_low]
target_psev_high  <- df_rect$psev_hi50[ df_rect$id == ind_psev_high]

s_sumrv_low  <- which.min(abs(scaled_sumrv[, ind_sumrv_low]  - target_sumrv_low))
s_sumrv_high <- which.min(abs(scaled_sumrv[, ind_sumrv_high] - target_sumrv_high))
s_psev_low   <- which.min(abs(posterior$psev[, ind_psev_low]  - target_psev_low))
s_psev_high  <- which.min(abs(posterior$psev[, ind_psev_high] - target_psev_high))

# Verify: print the target value and the actual sample value
cat("sumrv low  — ind:", ind_sumrv_low,
    " target Q25:", round(target_sumrv_low, 3),
    " sample value:", round(scaled_sumrv[s_sumrv_low, ind_sumrv_low], 3),
    " sample:", s_sumrv_low, "\n")

cat("sumrv high — ind:", ind_sumrv_high,
    " target Q75:", round(target_sumrv_high, 3),
    " sample value:", round(scaled_sumrv[s_sumrv_high, ind_sumrv_high], 3),
    " sample:", s_sumrv_high, "\n")

cat("psev low   — ind:", ind_psev_low,
    " target Q25:", round(target_psev_low, 4),
    " sample value:", round(posterior$psev[s_psev_low, ind_psev_low], 4),
    " sample:", s_psev_low, "\n")

cat("psev high  — ind:", ind_psev_high,
    " target Q75:", round(target_psev_high, 4),
    " sample value:", round(posterior$psev[s_psev_high, ind_psev_high], 4),
    " sample:", s_psev_high, "\n")


if (FALSE) {
s_sumrv_low  <- find_representative_sample(ind_sumrv_low,  "sumrv", sumrv_q25)
s_sumrv_high <- find_representative_sample(ind_sumrv_high, "sumrv", sumrv_q75)
s_psev_low   <- find_representative_sample(ind_psev_low,   "psev",  psev_q25)
s_psev_high  <- find_representative_sample(ind_psev_high,  "psev",  psev_q75)
}

extract_params <- function(s, i) {
  list(
    individual = i,
    sample     = s,
    beta       = beta_ind[s, i],
    gamma      = gamma_ind[s, i],
    omega      = omega_ind[s, i],
    pop        = pop_ind[s, i],
    Y0f        = Y0f[s],
    sumrv      = scaled_sumrv[s, i],
    psev       = posterior$psev[s, i]
  )
}

# ============================================================
# Package into params_list for run_sir_full
# ============================================================
params_list <- list(
  sumrv_low  = extract_params(s_sumrv_low,  ind_sumrv_low),
  sumrv_high = extract_params(s_sumrv_high, ind_sumrv_high),
  psev_low   = extract_params(s_psev_low,   ind_psev_low),
  psev_high  = extract_params(s_psev_high,  ind_psev_high)
)

# Quick sanity check
params_list %>%
  map_dfr(~ tibble(
    individual = .x$individual,
    sample     = .x$sample,
    beta       = round(.x$beta,  4),
    gamma      = round(.x$gamma, 4),
    omega      = round(.x$omega, 4),
    pop        = round(.x$pop,   0),
    sumrv      = round(.x$sumrv, 3),
    psev       = round(.x$psev,  5)
  ), .id = "scenario") %>%
  print()

params_list %>%
  map_dfr(~ tibble(
    individual = .x$individual,
    sample     = .x$sample,
    beta       = round(.x$beta,  4),
    gamma      = round(.x$gamma, 4),
    omega      = round(.x$omega, 4),
    R0         = round(.x$beta / (.x$gamma + .x$omega), 3),
    pop        = round(.x$pop,   0),
    sumrv      = round(.x$sumrv, 3),
    psev       = round(.x$psev,  5)
  ), .id = "scenario") %>%
  print()

##


run_sir_full <- function(p, tf = 63) {
  times  <- seq(0, tf, by = 1)
  state  <- c(S  = p$pop - p$Y0f,
              I  = p$Y0f,
              R  = 0,
              Cd = 0,
              Cc = 0)
  params <- c(beta  = p$beta,
              gamma = p$gamma,
              omega = p$omega,
              pop   = p$pop)
  out <- tryCatch(
    as.data.frame(ode(y = state, times = times,
                      func = sir_ode, parms = params,
                      method = "lsoda")),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  out %>% mutate(individual = p$individual,
                 sample     = p$sample,
                 sumrv      = p$sumrv,
                 psev       = p$psev)
}


df_sir_targets <- imap_dfr(params_list, function(p, scenario_name) {
  out <- run_sir_full(p, tf = 63)
  if (!is.null(out)) out %>% mutate(scenario = scenario_name)
})

scenario_labels <- c(
  sumrv_low  = glue::glue("Ind {ind_sumrv_low}  — sumrv Q25"),
  sumrv_high = glue::glue("Ind {ind_sumrv_high} — sumrv Q75"),
  psev_low   = glue::glue("Ind {ind_psev_low}   — psev Q25"),
  psev_high  = glue::glue("Ind {ind_psev_high}  — psev Q75")
)

scenario_colours <- c(
  sumrv_low  = "#92c5de",
  sumrv_high = "#2166ac",
  psev_low   = "#f4a582",
  psev_high  = "#b2182b"
)

df_sir_targets <- df_sir_targets %>%
  mutate(
    label    = scenario_labels[scenario],
    scenario = factor(scenario, levels = names(scenario_labels))
  )

df_pop_obs <- tibble(
  t      = mydata$tc_obs,
  cases  = cumsum(mydata$ncases),
  deaths = cumsum(mydata$ndeaths)
)

target_inds <- c(ind_sumrv_low, ind_sumrv_high,
                 ind_psev_low,  ind_psev_high)

ind_to_scenario <- c(
  setNames(rep("sumrv_low",  length(ind_sumrv_low)),  ind_sumrv_low),
  setNames(rep("sumrv_high", length(ind_sumrv_high)), ind_sumrv_high),
  setNames(rep("psev_low",   length(ind_psev_low)),   ind_psev_low),
  setNames(rep("psev_high",  length(ind_psev_high)),  ind_psev_high)
)

df_obs_targets <- tibble(
  individual = mydata$ind,
  t          = mydata$t,
  y_obs      = mydata$y_obs
) %>%
  filter(individual %in% target_inds) %>%
  mutate(scenario = ind_to_scenario[as.character(individual)])

df_vt_targets <- df_vt %>%
  filter(individual %in% target_inds) %>%
  mutate(scenario = ind_to_scenario[as.character(individual)])


df_sir_incident <- df_sir_targets %>%
  filter(time %in% mydata$tc_obs) %>%
  arrange(scenario, time) %>%
  group_by(scenario) %>%
  mutate(
    ncases  = c(first(Cc), diff(Cc)),   # incident cases per interval
    ndeaths = c(first(Cd), diff(Cd))    # incident deaths per interval
  ) %>%
  ungroup()

# Observed incident data
df_pop_incident <- tibble(
  t       = mydata$tc_obs,
  ncases  = mydata$ncases,
  ndeaths = mydata$ndeaths
)


n_samples <- length(posterior$beta_trans)

##
if (TRUE) {
df_ncases_all <- map_dfr(seq_len(n_ind), function(ind) {
  map_dfr(seq_len(n_samples), function(s) {
    state  <- c(S  = pop_ind[s, ind] - Y0f[s],
                I  = Y0f[s],
                R  = 0, Cd = 0, Cc = 0)
    params <- c(beta  = beta_ind[s, ind],
                gamma = gamma_ind[s, ind],
                omega = omega_ind[s, ind],
                pop   = pop_ind[s, ind])
    out <- tryCatch(
      as.data.frame(ode(y = state, times = c(0, mydata$tc_obs),
                        func = sir_ode, parms = params,
                        method = "lsoda")),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    # incident cases at each tc_obs
    tibble(
      sample   = s,
      time     = mydata$tc_obs,
      ncases   = diff(out$Cc)   # diff drops time=0
    )
  }) %>%
    group_by(time) %>%
    summarise(ncases_median = median(ncases, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(individual = ind)
}, .progress = TRUE)
}

obs_vec <- setNames(mydata$ncases, mydata$tc_obs)

df_distance <- df_ncases_all %>%
  group_by(individual) %>%
  summarise(
    # L2 distance (sum of squared differences)
    dist_l2    = sum((ncases_median - obs_vec[as.character(time)])^2,
                     na.rm = TRUE),
    # L1 distance (sum of absolute differences) — robust alternative
    dist_l1    = sum(abs(ncases_median - obs_vec[as.character(time)]),
                     na.rm = TRUE),
    # max absolute deviation
    dist_linf  = max(abs(ncases_median - obs_vec[as.character(time)]),
                     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dist_l2) %>%
  left_join(df_meta, by = c("individual" = "id"))

# Print top 5 closest by each metric
cat("=== Top 5 closest by L2 (sum of squared differences) ===\n")
df_distance %>%
  select(individual, sex_label, age_label, dist_l2, dist_l1, dist_linf) %>%
  slice_head(n = 5) %>%
  print()

cat("\n=== Top 5 closest by L1 (sum of absolute differences) ===\n")
df_distance %>%
  arrange(dist_l1) %>%
  select(individual, sex_label, age_label, dist_l2, dist_l1, dist_linf) %>%
  slice_head(n = 5) %>%
  print()

# Best individual by L2
ind_best <- df_distance %>% slice_min(dist_l2, n = 1) %>% pull(individual)
cat("\nBest individual (L2):", ind_best, "\n")

df_best <- df_ncases_all %>%
  filter(individual == ind_best)

p_check <- ggplot() +
  geom_line(
    data = df_best,
    aes(x = time, y = ncases_median),
    colour = "#2166ac", linewidth = 1.0
  ) +
  geom_line(
    data = tibble(t = mydata$tc_obs, ncases = mydata$ncases),
    aes(x = t, y = ncases),
    colour = "black", linewidth = 0.6, linetype = "dashed"
  ) +
  geom_point(
    data = tibble(t = mydata$tc_obs, ncases = mydata$ncases),
    aes(x = t, y = ncases),
    shape = 21, size = 2.2,
    fill = "white", colour = "black", stroke = 0.7
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x       = "Time (days)",
    y       = "Incident cases",
    title   = glue::glue("Best fitting individual: Ind {ind_best}"),
    subtitle = glue::glue(
      "L2 = {round(df_distance$dist_l2[df_distance$individual == ind_best], 1)}  |  ",
      "L1 = {round(df_distance$dist_l1[df_distance$individual == ind_best], 1)}  |  ",
      "{df_distance$sex_label[df_distance$individual == ind_best]}, ",
      "{df_distance$age_label[df_distance$individual == ind_best]}"
    )
  ) +
  my_theme_pub

p_check


###

best_colour  <- "#1a7a2e"   # green — distinct from existing palette
best_label   <- glue::glue("Ind {ind_best} — best fit")
best_meta    <- df_meta %>% filter(id == ind_best)
best_R0      <- beta_ind[, ind_best] / (gamma_ind[, ind_best] + omega_ind[, ind_best])
best_ann     <- glue::glue(
  "{best_meta$sex_label}, {best_meta$age_label}\n",
  "R\u2080 = {round(median(best_R0), 2)} ",
  "[{round(quantile(best_R0, 0.25), 2)}, {round(quantile(best_R0, 0.75), 2)}]"
)

if (FALSE) {
cat("Running", n_samples, "samples for best individual (ind", ind_best, ")...\n")
df_sir_best <- run_sir_individual(ind_best, "best") %>%
  mutate(scenario = "best")

df_sir_best_summary <- df_sir_best %>%
  filter(time %in% mydata$tc_obs) %>%
  arrange(sample, time) %>%
  group_by(sample) %>%
  mutate(
    ncases  = c(first(Cc), diff(Cc)),
    ndeaths = c(first(Cd), diff(Cd))
  ) %>%
  ungroup() %>%
  group_by(time) %>%
  summarise(
    ncases_median  = median(ncases,  na.rm = TRUE),
    ncases_lo50    = quantile(ncases,  0.25,  na.rm = TRUE),
    ncases_hi50    = quantile(ncases,  0.75,  na.rm = TRUE),
    ncases_lo95    = quantile(ncases,  0.025, na.rm = TRUE),
    ncases_hi95    = quantile(ncases,  0.975, na.rm = TRUE),
    ndeaths_median = median(ndeaths, na.rm = TRUE),
    ndeaths_lo50   = quantile(ndeaths, 0.25,  na.rm = TRUE),
    ndeaths_hi50   = quantile(ndeaths, 0.75,  na.rm = TRUE),
    ndeaths_lo95   = quantile(ndeaths, 0.025, na.rm = TRUE),
    ndeaths_hi95   = quantile(ndeaths, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(scenario = "best")

df_vt_best <- df_vt %>%
  filter(individual == ind_best, vt_median > 0) %>%
  mutate(
    scenario = "best",
    vt_lo50  = pmax(vt_lo50, 0),
    vt_hi50  = pmax(vt_hi50, 0),
    vt_lo95  = pmax(vt_lo95, 0),
    vt_hi95  = pmax(vt_hi95, 0)
  )

df_obs_best <- tibble(
  individual = mydata$ind,
  t          = mydata$t,
  y_obs      = mydata$y_obs
) %>%
  filter(individual == ind_best, y_obs > 0) %>%
  mutate(scenario = "best")

}

# Figure 3 here

library(glue)
library(cowplot)


my_theme_pub <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.text        = element_text(colour = "black", size = 9),
    plot.title       = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    legend.position  = "bottom"
  )

# Per-individual posterior medians
df_ind_medians <- tibble(
  id           = seq_len(n_ind),
  psev_median  = apply(posterior$psev, 2, median),
  sumrv_median = apply(scaled_sumrv,   2, median)
)

ind_psev_low   <- df_ind_medians %>% slice_min(psev_median,  n = 1) %>% pull(id)
ind_psev_high  <- df_ind_medians %>% slice_max(psev_median,  n = 1) %>% pull(id)
ind_sumrv_low  <- df_ind_medians %>% slice_min(sumrv_median, n = 1) %>% pull(id)
ind_sumrv_high <- df_ind_medians %>% slice_max(sumrv_median, n = 1) %>% pull(id)

target_inds <- c(
  psev_low   = ind_psev_low,
  psev_high  = ind_psev_high,
  sumrv_low  = ind_sumrv_low,
  sumrv_high = ind_sumrv_high
)

cat("Target individuals:\n")
cat("  psev  lowest  — Ind", ind_psev_low,
    " median:", round(df_ind_medians$psev_median[ind_psev_low],   5), "\n")
cat("  psev  highest — Ind", ind_psev_high,
    " median:", round(df_ind_medians$psev_median[ind_psev_high],  5), "\n")
cat("  sumrv lowest  — Ind", ind_sumrv_low,
    " median:", round(df_ind_medians$sumrv_median[ind_sumrv_low],  3), "\n")
cat("  sumrv highest — Ind", ind_sumrv_high,
    " median:", round(df_ind_medians$sumrv_median[ind_sumrv_high], 3), "\n")


obs_vec <- setNames(mydata$ncases, mydata$tc_obs)

best_colour <- "#1a7a2e"

sc_levels_5 <- c(names(target_inds), "best")

scenario_colours_5 <- c(
  psev_low   = "#92c5de",
  psev_high  = "#2166ac",
  sumrv_low  = "#f4a582",
  sumrv_high = "#b2182b",
  best       = best_colour
)

# Strip labels: Ind # only
scenario_labels_5 <- c(
  psev_low   = glue("Ind {target_inds['psev_low']}"),
  psev_high  = glue("Ind {target_inds['psev_high']}"),
  sumrv_low  = glue("Ind {target_inds['sumrv_low']}"),
  sumrv_high = glue("Ind {target_inds['sumrv_high']}"),
  best       = glue("Ind {ind_best}")
)

sc_labeller_5 <- labeller(scenario = scenario_labels_5)

# In-panel annotations: sex + age only
df_annotations_5 <- bind_rows(
  map_dfr(names(target_inds), function(sc) {
    meta <- df_meta %>% filter(id == target_inds[sc])
    tibble(scenario = sc,
           ann_text = glue("{meta$sex_label}, {meta$age_label}"))
  }),
  tibble(scenario = "best",
         ann_text = glue("{best_meta$sex_label}, {best_meta$age_label}"))
) %>%
  mutate(scenario = factor(scenario, levels = sc_levels_5))

run_sir_individual <- function(ind, scenario_name) {
  map_dfr(seq_len(n_samples), function(s) {
    state  <- c(S = pop_ind[s, ind] - Y0f[s],
                I = Y0f[s], R = 0, Cd = 0, Cc = 0)
    params <- c(beta  = beta_ind[s, ind],
                gamma = gamma_ind[s, ind],
                omega = omega_ind[s, ind],
                pop   = pop_ind[s, ind])
    out <- tryCatch(
      as.data.frame(ode(y = state, times = seq(0, 63, by = 1),
                        func = sir_ode, parms = params, method = "lsoda")),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    out %>% mutate(sample = s)
  })
}

summarise_sir <- function(df_raw, scenario_name) {
  df_raw %>%
    filter(time %in% mydata$tc_obs) %>%
    arrange(sample, time) %>%
    group_by(sample) %>%
    mutate(ncases  = c(first(Cc), diff(Cc)),
           ndeaths = c(first(Cd), diff(Cd))) %>%
    ungroup() %>%
    group_by(time) %>%
    summarise(
      ncases_median  = median(ncases,  na.rm = TRUE),
      ncases_lo50    = quantile(ncases,  0.25,  na.rm = TRUE),
      ncases_hi50    = quantile(ncases,  0.75,  na.rm = TRUE),
      ncases_lo95    = quantile(ncases,  0.025, na.rm = TRUE),
      ncases_hi95    = quantile(ncases,  0.975, na.rm = TRUE),
      ndeaths_median = median(ndeaths, na.rm = TRUE),
      ndeaths_lo50   = quantile(ndeaths, 0.25,  na.rm = TRUE),
      ndeaths_hi50   = quantile(ndeaths, 0.75,  na.rm = TRUE),
      ndeaths_lo95   = quantile(ndeaths, 0.025, na.rm = TRUE),
      ndeaths_hi95   = quantile(ndeaths, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_name)
}

all_inds_5 <- c(target_inds, best = ind_best)

df_sir_all5 <- imap_dfr(all_inds_5, function(ind, sc) {
  cat("Running", n_samples, "samples for", sc, "(Ind", ind, ")...\n")
  run_sir_individual(ind, sc) %>%
    summarise_sir(sc)
}) %>%
  mutate(scenario = factor(scenario, levels = sc_levels_5))

#####

library(glue)
library(cowplot)

my_theme_pub <- theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.text        = element_text(colour = "black", size = 9),
    plot.title       = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    legend.position  = "bottom"
  )

# Per-individual posterior medians
df_ind_medians <- tibble(
  id           = seq_len(n_ind),
  psev_median  = apply(posterior$psev, 2, median),
  sumrv_median = apply(scaled_sumrv,   2, median)
)

ind_psev_low   <- df_ind_medians %>% slice_min(psev_median,  n = 1) %>% pull(id)
ind_psev_high  <- df_ind_medians %>% slice_max(psev_median,  n = 1) %>% pull(id)
ind_sumrv_low  <- df_ind_medians %>% slice_min(sumrv_median, n = 1) %>% pull(id)
ind_sumrv_high <- df_ind_medians %>% slice_max(sumrv_median, n = 1) %>% pull(id)

target_inds <- c(
  psev_low   = ind_psev_low,
  psev_high  = ind_psev_high,
  sumrv_low  = ind_sumrv_low,
  sumrv_high = ind_sumrv_high
)

cat("Target individuals:\n")
cat("  psev  lowest  — Ind", ind_psev_low,
    " median:", round(df_ind_medians$psev_median[ind_psev_low],   5), "\n")
cat("  psev  highest — Ind", ind_psev_high,
    " median:", round(df_ind_medians$psev_median[ind_psev_high],  5), "\n")
cat("  sumrv lowest  — Ind", ind_sumrv_low,
    " median:", round(df_ind_medians$sumrv_median[ind_sumrv_low],  3), "\n")
cat("  sumrv highest — Ind", ind_sumrv_high,
    " median:", round(df_ind_medians$sumrv_median[ind_sumrv_high], 3), "\n")

# Run ODE for all individuals at tc_obs time points only

# use mean of pop_ind

mypop_ind = mean(pop_ind)
myinitialI0 = mean(Y0f)

df_ncases_all <- map_dfr(seq_len(n_ind), function(ind) {
  map_dfr(seq_len(n_samples), function(s) {
    state  <- c(S = mypop_ind - myinitialI0,
                I = myinitialI0, R = 0, Cd = 0, Cc = myinitialI0)
    params <- c(beta  = beta_ind[s, ind],
                gamma = gamma_ind[s, ind],
                omega = omega_ind[s, ind],
                pop   = pop_ind[s, ind])
    out <- tryCatch(
      as.data.frame(ode(y = state, times = c(0, mydata$tc_obs),
                        func = sir_ode, parms = params, method = "lsoda")),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    tibble(sample = s, time = mydata$tc_obs, ncases = diff(out$Cc))
  }) %>%
    group_by(time) %>%
    summarise(ncases_median = median(ncases, na.rm = TRUE), .groups = "drop") %>%
    mutate(individual = ind)
}, .progress = TRUE)

obs_vec <- setNames(mydata$ncases, mydata$tc_obs)

df_distance <- df_ncases_all %>%
  group_by(individual) %>%
  summarise(
    dist_l2   = sum((ncases_median - obs_vec[as.character(time)])^2,  na.rm = TRUE),
    dist_l1   = sum(abs(ncases_median - obs_vec[as.character(time)]), na.rm = TRUE),
    dist_linf = max(abs(ncases_median - obs_vec[as.character(time)]), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dist_l2) %>%
  left_join(df_meta, by = c("individual" = "id"))

ind_best  <- df_distance %>% slice_min(dist_l2, n = 1) %>% pull(individual)
best_meta <- df_meta %>% filter(id == ind_best)
cat("Best-fitting individual (L2): Ind", ind_best, "\n")

best_colour <- "#1a7a2e"

sc_levels_5 <- c(names(target_inds), "best")

scenario_colours_5 <- c(
  psev_low   = "#92c5de",
  psev_high  = "#2166ac",
  sumrv_low  = "#f4a582",
  sumrv_high = "#b2182b",
  best       = best_colour
)

# Strip labels: Ind # only
scenario_labels_5 <- c(
  psev_low   = glue("Ind {target_inds['psev_low']}"),
  psev_high  = glue("Ind {target_inds['psev_high']}"),
  sumrv_low  = glue("Ind {target_inds['sumrv_low']}"),
  sumrv_high = glue("Ind {target_inds['sumrv_high']}"),
  best       = glue("Ind {ind_best}")
)

sc_labeller_5 <- labeller(scenario = scenario_labels_5)

# In-panel annotations: sex + age only
df_annotations_5 <- bind_rows(
  map_dfr(names(target_inds), function(sc) {
    meta <- df_meta %>% filter(id == target_inds[sc])
    tibble(scenario = sc,
           ann_text = glue("{meta$sex_label}, {meta$age_label}"))
  }),
  tibble(scenario = "best",
         ann_text = glue("{best_meta$sex_label}, {best_meta$age_label}"))
) %>%
  mutate(scenario = factor(scenario, levels = sc_levels_5))

run_sir_individual <- function(ind, scenario_name) {
  map_dfr(seq_len(n_samples), function(s) {
    state  <- c(S = pop_ind[s, ind] - Y0f[s],
                I = Y0f[s], R = 0, Cd = 0, Cc = myinitialI0)
    params <- c(beta  = beta_ind[s, ind],
                gamma = gamma_ind[s, ind],
                omega = omega_ind[s, ind],
                pop   = pop_ind[s, ind])
    out <- tryCatch(
      as.data.frame(ode(y = state, times = seq(0, 63, by = 1),
                        func = sir_ode, parms = params, method = "lsoda")),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    out %>% mutate(sample = s)
  })
}

summarise_sir <- function(df_raw, scenario_name) {
  df_raw %>%
    filter(time %in% mydata$tc_obs) %>%
    arrange(sample, time) %>%
    group_by(sample) %>%
    mutate(ncases  = c(first(Cc), diff(Cc)),
           ndeaths = c(first(Cd), diff(Cd))) %>%
    ungroup() %>%
    group_by(time) %>%
    summarise(
      ncases_median  = median(ncases,  na.rm = TRUE),
      ncases_lo50    = quantile(ncases,  0.25,  na.rm = TRUE),
      ncases_hi50    = quantile(ncases,  0.75,  na.rm = TRUE),
      ncases_lo95    = quantile(ncases,  0.025, na.rm = TRUE),
      ncases_hi95    = quantile(ncases,  0.975, na.rm = TRUE),
      ndeaths_median = median(ndeaths, na.rm = TRUE),
      ndeaths_lo50   = quantile(ndeaths, 0.25,  na.rm = TRUE),
      ndeaths_hi50   = quantile(ndeaths, 0.75,  na.rm = TRUE),
      ndeaths_lo95   = quantile(ndeaths, 0.025, na.rm = TRUE),
      ndeaths_hi95   = quantile(ndeaths, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(scenario = scenario_name)
}

all_inds_5 <- c(target_inds, best = ind_best)
all_inds_5 <- c(target_inds)

sc_levels <- names(target_inds)
sc_levels_5 <- c(names(target_inds), "best")

df_sir_all5 <- imap_dfr(all_inds_5, function(ind, sc) {
  cat("Running", n_samples, "samples for", sc, "(Ind", ind, ")...\n")
  run_sir_individual(ind, sc) %>%
    summarise_sir(sc)
}) %>%
  mutate(scenario = factor(scenario, levels = sc_levels))

df_sir_pop_raw <- map_dfr(seq_len(n_samples), function(s) {
  p      <- list(beta  = posterior$beta_trans[s],
                 gamma = posterior$gamma_trans[s],
                 omega = posterior$omega_trans[s],
                 pop   = mean(pop_ind[s, ]),
                 Y0f   = Y0f[s])
  state  <- c(S = p$pop - p$Y0f, I = p$Y0f, R = 0, Cd = 0, Cc = myinitialI0)
  params <- c(beta = p$beta, gamma = p$gamma,
              omega = p$omega, pop = p$pop)
  out <- tryCatch(
    as.data.frame(ode(y = state, times = c(0, mydata$tc_obs),
                      func = sir_ode, parms = params, method = "lsoda")),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  tibble(sample = s, time = mydata$tc_obs,
         ncases  = diff(out$Cc),
         ndeaths = diff(out$Cd))
}, .progress = TRUE)

df_pop_post <- df_sir_pop_raw %>%
  group_by(time) %>%
  summarise(
    ncases_median  = median(ncases,  na.rm = TRUE),
    ncases_lo50    = quantile(ncases,  0.25,  na.rm = TRUE),
    ncases_hi50    = quantile(ncases,  0.75,  na.rm = TRUE),
    ncases_lo95    = quantile(ncases,  0.025, na.rm = TRUE),
    ncases_hi95    = quantile(ncases,  0.975, na.rm = TRUE),
    ndeaths_median = median(ndeaths, na.rm = TRUE),
    ndeaths_lo50   = quantile(ndeaths, 0.25,  na.rm = TRUE),
    ndeaths_hi50   = quantile(ndeaths, 0.75,  na.rm = TRUE),
    ndeaths_lo95   = quantile(ndeaths, 0.025, na.rm = TRUE),
    ndeaths_hi95   = quantile(ndeaths, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

df_pop_incident <- tibble(
  t       = mydata$tc_obs,
  ncases  = mydata$ncases,
  ndeaths = mydata$ndeaths
)

sc_levels <- names(target_inds)
sc_levels_5 <- c(names(target_inds), "best")

df_vt_all5 <- imap_dfr(all_inds_5, function(ind, sc) {
  df_vt %>%
    filter(individual == ind, vt_median > 0) %>%
    mutate(
      scenario = sc,
      vt_lo50  = pmax(vt_lo50, 0),
      vt_hi50  = pmax(vt_hi50, 0),
      vt_lo95  = pmax(vt_lo95, 0),
      vt_hi95  = pmax(vt_hi95, 0)
    )
}) %>%
  mutate(scenario = factor(scenario, levels = sc_levels_5))

df_obs_all5 <- imap_dfr(all_inds_5, function(ind, sc) {
  tibble(
    individual = mydata$ind,
    t          = mydata$t,
    y_obs      = mydata$y_obs
  ) %>%
    filter(individual == ind, y_obs > 0) %>%
    mutate(scenario = sc)
}) %>%
  mutate(scenario = factor(scenario, levels = sc_levels_5))


df_vt_all5_old <- df_vt_all5
df_vt_all5_old %>%
  filter(scenario!= "best") -> df_vt_all5

df_obs_all5_old <- df_obs_all5
df_obs_all5_old %>%
  filter(scenario!= "best") -> df_obs_all5

scenario_colours_5
scenario_colours

df_annotations_5_old <- df_annotations_5
df_annotations_5_old %>%
  filter(scenario!="best") -> df_annotations_5

my_theme_pub_y <- theme_classic(base_size = 16, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 16),
    axis.title       = element_text(size = 15),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 16),
    plot.margin      = margin(8, 8, 8, 8)
  )


# esse que vale
# ── Row 1: Kinetics ─────────────────────────────────────────
p_row_kinetics <- ggplot() +
  geom_ribbon(
    data = df_vt_all5,
    aes(x = t, ymin = vt_lo95/log(10), ymax = vt_hi95/log(10), fill = scenario),
    alpha = 0.12
  ) +
  geom_ribbon(
    data = df_vt_all5,
    aes(x = t, ymin = vt_lo50/log(10), ymax = vt_hi50/log(10), fill = scenario),
    alpha = 0.30
  ) +
  geom_line(
    data = df_vt_all5,
    aes(x = t, y = vt_median/log(10), colour = scenario),
    linewidth = 0.9
  ) +
  geom_point(
    data = df_obs_all5,
    aes(x = t, y = y_obs/log(10), colour = scenario),
    shape = 21, size = 1.8, fill = "white", stroke = 0.7
  ) +
#  geom_text(
#    data = df_annotations_5,
#    aes(label = ann_text, colour = scenario),
#    x = -Inf, y = Inf,
#    hjust = -0.08, vjust = 1.3,
#    size = 2.6, lineheight = 1.2,
#    show.legend = FALSE
#  ) +
  scale_colour_manual(values = scenario_colours,
                      labels = scenario_labels, name = NULL) +
  scale_fill_manual(  values = scenario_colours, guide = "none") +
#  facet_wrap(~ scenario, nrow = 1, labeller = sc_labeller) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(x = NULL, y = "Viral titer", title = "Viral kinetics",
       tag = "A") +
  my_theme_pub_y +
  theme(legend.position  = "none",
#        strip.text       = element_text(size = 9, face = "bold"),
        strip.text       = element_blank(),
        strip.background = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

p_row_kinetics

cases_max  <- max(df_sir_all5$ncases_hi95,
                  df_pop_post$ncases_hi95,
                  mydata$ncases,
                  na.rm = TRUE)

deaths_max <- max(df_sir_all5$ndeaths_hi95,
                  df_pop_post$ndeaths_hi95,
                  mydata$ndeaths,
                  na.rm = TRUE)

cases_lim  <- c(0, cases_max  * 1.05)
deaths_lim <- c(0, deaths_max * 1.05)

p_row_cases <- ggplot() +
  geom_ribbon(
    data = df_sir_all5,
    aes(x = time, ymin = ncases_lo95, ymax = ncases_hi95, fill = scenario),
    alpha = 0.12
  ) +
  geom_ribbon(
    data = df_sir_all5,
    aes(x = time, ymin = ncases_lo50, ymax = ncases_hi50, fill = scenario),
    alpha = 0.30
  ) +
  geom_line(
    data = df_sir_all5,
    aes(x = time, y = ncases_median, colour = scenario),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = scenario_colours_5,
                      labels = scenario_labels_5, name = NULL) +
  scale_fill_manual(  values = scenario_colours_5, guide = "none") +
  scale_y_continuous(labels = scales::comma, limits = cases_lim) +
  facet_wrap(~ scenario, nrow = 1, labeller = sc_labeller_5) +
  labs(x = NULL, y = "Incident cases", title = "Cases per interval",
       tag = "B") +
  my_theme_pub_y +
  theme(legend.position  = "none",
        strip.text       = element_blank(),
        axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank())

p_row_cases

#df_sir_all5_old <- df_sir_all5
#df_annotations_5_old %>%
#  filter(scenario!="best") -> df_sir_all5

# Strip labels: Ind # only
if (FALSE) {
  scenario_labels_5 <- c(
  psev_low   = glue("Ind {target_inds['psev_low']}"),
  psev_high  = glue("Ind {target_inds['psev_high']}"),
  sumrv_low  = glue("Ind {target_inds['sumrv_low']}"),
  sumrv_high = glue("Ind {target_inds['sumrv_high']}"),
  best       = glue("Ind {ind_best}")
)
}
  
#sc_labeller_5 <- labeller(scenario = scenario_labels_5)
sc_labeller <- labeller(scenario = scenario_labels)


# ── Row 3: Deaths ────────────────────────────────────────────
p_row_deaths <- ggplot() +
  geom_ribbon(
    data = df_sir_all5,
    aes(x = time, ymin = ndeaths_lo95, ymax = ndeaths_hi95, fill = scenario),
    alpha = 0.12
  ) +
  geom_ribbon(
    data = df_sir_all5,
    aes(x = time, ymin = ndeaths_lo50, ymax = ndeaths_hi50, fill = scenario),
    alpha = 0.30
  ) +
  geom_line(
    data = df_sir_all5,
    aes(x = time, y = ndeaths_median, colour = scenario),
    linewidth = 0.9
  ) +
  scale_colour_manual(values = scenario_colours,
                      labels = scenario_labels, name = NULL) +
  scale_fill_manual(  values = scenario_colours, guide = "none") +
  scale_y_continuous(labels = scales::comma, limits = deaths_lim) +
  facet_wrap(~ scenario, nrow = 1, labeller = sc_labeller) +
  labs(x = "Time (days)", y = "Incident deaths",
       title = "Deaths per interval",
       tag   = "C") +
  my_theme_pub_y +
  theme(legend.position = "none",
        strip.text      = element_blank())

p_row_deaths

# ── Population side column ───────────────────────────────────
p_pop_cases <- ggplot() +
  geom_ribbon(
    data = df_pop_post,
    aes(x = time, ymin = ncases_lo95, ymax = ncases_hi95),
    fill = "#cccccc", alpha = 0.35
  ) +
  geom_ribbon(
    data = df_pop_post,
    aes(x = time, ymin = ncases_lo50, ymax = ncases_hi50),
    fill = "#888888", alpha = 0.35
  ) +
  geom_line(
    data = df_pop_post,
    aes(x = time, y = ncases_median),
    colour = "#333333", linewidth = 0.9
  ) +
  geom_line(
    data = df_pop_incident,
    aes(x = t, y = ncases),
    colour = "black", linewidth = 0.6, linetype = "dashed"
  ) +
  scale_y_continuous(labels = scales::comma, limits = cases_lim) +
  labs(x = NULL, y = "Incident cases",
       tag   = "D",
       title = "Population level\n(\u03b2/\u03b3/\u03c9_trans)") +
  my_theme_pub_y +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title   = element_text(size = 9, face = "bold", hjust = 0.5))

p_pop_deaths <- ggplot() +
  geom_ribbon(
    data = df_pop_post,
    aes(x = time, ymin = ndeaths_lo95, ymax = ndeaths_hi95),
    fill = "#cccccc", alpha = 0.35
  ) +
  geom_ribbon(
    data = df_pop_post,
    aes(x = time, ymin = ndeaths_lo50, ymax = ndeaths_hi50),
    fill = "#888888", alpha = 0.35
  ) +
  geom_line(
    data = df_pop_post,
    aes(x = time, y = ndeaths_median),
    colour = "#333333", linewidth = 0.9
  ) +
  geom_line(
    data = df_pop_incident,
    aes(x = t, y = ndeaths),
    colour = "black", linewidth = 0.6, linetype = "dashed"
  ) +
  scale_y_continuous(labels = scales::comma, limits = deaths_lim) +
  labs(x = "Time (days)", y = "Incident deaths", title = NULL) +
  my_theme_pub_y

p_pop_col <- wrap_elements(
  full = (p_pop_cases / p_pop_deaths) + plot_layout(heights = c(1, 1))
) +
  theme(plot.background = element_rect(
    fill = "#f7f7f7", colour = "#bbbbbb", linewidth = 0.6
  ))

# ── Shared legend ─────────────────────────────────────────────
leg <- get_legend(
  ggplot(df_sir_all5,
         aes(x = time, y = ncases_median, colour = scenario)) +
    geom_line() +
    scale_colour_manual(values = scenario_colours_5,
                        labels = scenario_labels_5, name = NULL) +
    my_theme_pub_y +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 14),
          legend.key.width = unit(1.2, "cm"))
)

p_pop_col <- (plot_spacer() / p_pop_cases / p_pop_deaths / plot_spacer()) +
  plot_layout(heights = c(1.2, 1, 1, 0.12)) +
  plot_annotation(tag_levels = list(c("", "D", "", ""))) +
  theme(plot.background = element_rect(
    fill = "#f7f7f7", colour = "#bbbbbb", linewidth = 0.6
  ))

p_pop_col <- (p_pop_cases / p_pop_deaths / plot_spacer()) +
  plot_layout(heights = c(2.2, 1, 0.12)) +  # 2.2 = kinetics (1.2) + cases (1)
  plot_annotation(tag_levels = list(c("D", "", ""))) +
  theme(plot.background = element_rect(
    fill = "#f7f7f7", colour = "#bbbbbb", linewidth = 0.6
  ))


p_pop_col <- (p_pop_cases / p_pop_deaths / plot_spacer()) +
  plot_layout(heights = c(1.6, 1.6, 0.12)) +
  plot_annotation(tag_levels = list(c("D", "", ""))) +
  theme(plot.background = element_rect(
    fill = "#f7f7f7", colour = "#bbbbbb", linewidth = 0.6
  ))

# No theme on p_pop_col itself — spacer stays transparent
p_pop_col <- (p_pop_cases / p_pop_deaths / plot_spacer()) +
  plot_layout(heights = c(1.6, 1.6, 0.12)) +
  plot_annotation(tag_levels = list(c("D", "", "")))

p_main_col <- (p_row_kinetics /
                 p_row_cases    /
                 p_row_deaths   /
                 plot_grid(leg)) +
  plot_layout(heights = c(1.2, 1, 1, 0.12))

p_final <- (p_main_col | p_pop_col) +
  plot_layout(widths = c(4, 1)) +
  plot_annotation(
#    caption = paste0(
#      "Cols 1–2: extreme psev (lowest/highest). ",
#      "Cols 3–4: extreme sumrv (lowest/highest). ",
#      "Col 5 (green): best-fitting individual to observed cases (L2).\n",
#      "Right panel (grey): population-level \u03b2/\u03b3/\u03c9_trans posterior. ",
#      "Dashed: observed data. Ribbons: 50% (dark) and 95% (light) credible intervals."
#    ),
    theme = theme(
   #   plot.caption      = element_text(size = 8, colour = "grey50", hjust = 0),
      plot.tag          = element_text(size = 13, face = "bold", family = "serif"),
      plot.tag.position = "topleft"
    )
  )

p_final

ggsave("sir_extreme_individuals.pdf",  plot = p_final,
       width = 20, height = 12, device = cairo_pdf)
ggsave("sir_extreme_individuals.tiff", plot = p_final,
       width = 20, height = 12, dpi = 600, compression = "lzw")

cat("Saved sir_extreme_individuals.pdf and .tiff\n")


fmt_miqr95 <- function(x, digits = 3) {
  glue::glue(
    "{round(median(x), digits)} ",
    "[{round(quantile(x, 0.25), digits)}, {round(quantile(x, 0.75), digits)}] ",
    "({round(quantile(x, 0.025), digits)}, {round(quantile(x, 0.975), digits)})"
  )
}

all_inds_table <- c(target_inds, best = ind_best)

df_target_table <- map_dfr(names(all_inds_table), function(sc) {
  ind  <- all_inds_table[sc]
  meta <- df_meta %>% filter(id == ind)
  R0   <- beta_ind[, ind] / (gamma_ind[, ind] + omega_ind[, ind])
  
  tibble(
    scenario  = sc,
    individual = as.character(ind),
    sex        = meta$sex_label,
    age        = meta$age_label,
    n_obs      = sum(mydata$ind == ind),
    beta_med   = round(median(beta_ind[,  ind]), 4),
    gamma_med  = round(median(gamma_ind[, ind]), 4),
    omega_med  = round(median(omega_ind[, ind]), 4),
    R0         = fmt_miqr95(R0,                    digits = 2),
    psev       = fmt_miqr95(posterior$psev[, ind], digits = 4)
  )
})

df_sir_totals <- imap_dfr(all_inds_table, function(ind, sc) {
  map_dfr(seq_len(n_samples), function(s) {
    state  <- c(S = pop_ind[s, ind] - Y0f[s],
                I = Y0f[s], R = 0, Cd = 0, Cc = myinitialI0)
    params <- c(beta  = beta_ind[s, ind],
                gamma = gamma_ind[s, ind],
                omega = omega_ind[s, ind],
                pop   = pop_ind[s, ind])
    out <- tryCatch(
      as.data.frame(ode(y = state, times = c(0, 63),
                        func = sir_ode, parms = params,
                        method = "lsoda")),
      error = function(e) NULL
    )
    if (is.null(out)) return(NULL)
    tibble(sample      = s,
           total_cases  = tail(out$Cc, 1),
           total_deaths = tail(out$Cd, 1))
  }) %>%
    summarise(
      cases_summary  = fmt_miqr95(total_cases,  digits = 0),
      deaths_summary = fmt_miqr95(total_deaths, digits = 0),
      .groups = "drop"
    ) %>%
    mutate(scenario = sc)
})

df_target_table <- df_target_table %>%
  left_join(df_sir_totals, by = "scenario") %>%
  mutate(
    scenario = factor(scenario,
                      levels = c(names(target_inds), "best"))
  ) %>%
  arrange(scenario)

row_colours <- c(
  psev_low   = "#92c5de",
  psev_high  = "#2166ac",
  sumrv_low  = "#f4a582",
  sumrv_high = "#b2182b",
  best       = "#1a7a2e"
)

row_text_colours <- c(
  psev_low   = "black",
  psev_high  = "white",
  sumrv_low  = "black",
  sumrv_high = "white",
  best       = "white"
)

tbl_targets <- df_target_table %>%
  select(scenario, individual, sex, age, n_obs,
         beta_med, gamma_med, omega_med,
         R0, psev, cases_summary, deaths_summary) %>%
  gt(rowname_col = "scenario") %>%
  # ── column labels ─────────────────────────────────────────
  cols_label(
    individual     = "Ind",
    sex            = "Sex",
    age            = "Age",
    n_obs          = "N obs",
    beta_med       = md("β_ind<br>median"),
    gamma_med      = md("γ_ind<br>median"),
    omega_med      = md("ω_ind<br>median"),
    R0             = md("R₀<br>median [IQR] (95%)"),
    psev           = md("*p*_sev<br>median [IQR] (95%)"),
    cases_summary  = md("Total cases<br>median [IQR] (95%)"),
    deaths_summary = md("Total deaths<br>median [IQR] (95%)")
  ) %>%
  # ── spanners ──────────────────────────────────────────────
  tab_spanner(
    label   = "Individual",
    columns = c(individual, sex, age, n_obs)
  ) %>%
  tab_spanner(
    label   = md("SIR parameters — posterior median"),
    columns = c(beta_med, gamma_med, omega_med)
  ) %>%
  tab_spanner(
    label   = md("Posterior summaries"),
    columns = c(R0, psev)
  ) %>%
  tab_spanner(
    label   = md("SIR outcomes (t = 63 days)"),
    columns = c(cases_summary, deaths_summary)
  ) %>%
  # ── header ────────────────────────────────────────────────
  tab_header(
    title    = md("**Target individuals — summary statistics**"),
    subtitle = md("Posterior summaries: median [IQR] (2.5%, 97.5%)")
  ) %>%
  # ── row stub colours ──────────────────────────────────────
  tab_style(
    style     = list(cell_fill(color = "#92c5de"),
                     cell_text(weight = "bold", color = "black")),
    locations = cells_stub(rows = scenario == "psev_low")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#2166ac"),
                     cell_text(weight = "bold", color = "white")),
    locations = cells_stub(rows = scenario == "psev_high")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#f4a582"),
                     cell_text(weight = "bold", color = "black")),
    locations = cells_stub(rows = scenario == "sumrv_low")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#b2182b"),
                     cell_text(weight = "bold", color = "white")),
    locations = cells_stub(rows = scenario == "sumrv_high")
  ) %>%
  tab_style(
    style     = list(cell_fill(color = "#1a7a2e"),
                     cell_text(weight = "bold", color = "white")),
    locations = cells_stub(rows = scenario == "best")
  ) %>%
  # ── column tints ──────────────────────────────────────────
  tab_style(
    style     = cell_fill(color = "#f7fbff"),
    locations = cells_body(columns = c(beta_med, gamma_med, omega_med))
  ) %>%
  tab_style(
    style     = cell_fill(color = "#fff7f0"),
    locations = cells_body(columns = c(R0, psev))
  ) %>%
  tab_style(
    style     = cell_fill(color = "#f4faf4"),
    locations = cells_body(columns = c(cases_summary, deaths_summary))
  ) %>%
  # ── alignment ─────────────────────────────────────────────
  cols_align(align = "center",
             columns = c(individual, sex, age, n_obs,
                         beta_med, gamma_med, omega_med)) %>%
  cols_align(align = "left",
             columns = c(R0, psev, cases_summary, deaths_summary)) %>%
  # ── footnotes ─────────────────────────────────────────────
  tab_footnote(
    footnote  = "β_ind, γ_ind, ω_ind: SIR transmission, recovery and mortality rates.",
    locations = cells_column_spanners("SIR parameters — posterior median")
  ) %>%
  tab_footnote(
    footnote  = "R₀ = β_ind / (γ_ind + ω_ind) per posterior sample.",
    locations = cells_column_labels(columns = R0)
  ) %>%
  tab_footnote(
    footnote  = "Best-fitting individual selected by minimum L2 distance to observed incident cases.",
    locations = cells_stub(rows = scenario == "best")
  ) %>%
  tab_options(
    table.font.names                  = "serif",
    table.font.size                   = 12,
    heading.align                     = "left",
    column_labels.border.top.width    = px(2),
    column_labels.border.bottom.width = px(2),
    table.border.top.width            = px(2),
    table.border.bottom.width         = px(2),
    data_row.padding                  = px(6),
    stub.border.width                 = px(0)
  )

gtsave(tbl_targets, "target_individuals_summary.html")
gtsave(tbl_targets, "target_individuals_summary.pdf")
gtsave(tbl_targets, "target_individuals_summary.rtf")

print(tbl_targets)

###

summarise_param <- function(mat, param_name) {
  apply(mat, 2, function(x) {
    c(median = median(x), q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
  }) |>
    t() |>
    as.data.frame() |>
    setNames(paste0(param_name, c("_median", "_q25", "_q75"))) |>
    mutate(id = 1:nrow(.))
}

summarise_param <- function(mat, param_name) {
  result <- apply(mat, 2, function(x) {
    c(median = median(x), q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
  }) |>
    t() |>
    as.data.frame() |>
    setNames(paste0(param_name, c("_median", "_q25", "_q75")))
  
  result$id <- seq_len(nrow(result))
  result
}

tbl <- left_join(
  summarise_param(scaled_sumrv, "sumrv"),
  summarise_param(posterior$psev,  "psev"),
  by = "id"
)

tbl |>
  mutate(
    sumrv_iqr = paste0(round(sumrv_median, 2), " [", round(sumrv_q25, 2), ", ", round(sumrv_q75, 2), "]"),
    psev_iqr  = paste0(round(psev_median,  2), " [", round(psev_q25,  2), ", ", round(psev_q75,  2), "]")
  ) |>
  select(id, sumrv_iqr, psev_iqr) |>
  gt() |>
  cols_label(
    id        = "ID",
    sumrv_iqr = "sumrv — Median [IQR]",
    psev_iqr  = "psev — Median [IQR]"
  ) |>
  tab_header(title = "Posterior summaries: sumrv and psev by individual")


library(gt)

tbl_gt <- tbl |>
  mutate(
    sumrv_iqr = paste0(round(sumrv_median, 2), " [", round(sumrv_q25, 2), ", ", round(sumrv_q75, 2), "]"),
    psev_iqr  = paste0(round(psev_median,  2), " [", round(psev_q25,  2), ", ", round(psev_q75,  2), "]")
  ) |>
  select(id, sumrv_iqr, psev_iqr) |>
  gt() |>
  cols_label(
    id        = "ID",
    sumrv_iqr = "sumrv — Median [IQR]",
    psev_iqr  = "psev — Median [IQR]"
  ) |>
  tab_header(title = "Posterior summaries: sumrv and psev by individual")

tbl_gt <- tbl |>
  left_join(df_distance |> select(individual, dist_l2), by = c("id"="individual")) |>
  mutate(
    sumrv_iqr = paste0(round(sumrv_median, 2), " [", round(sumrv_q25, 2), ", ", round(sumrv_q75, 2), "]"),
    psev_iqr  = paste0(round(psev_median,  2), " [", round(psev_q25,  2), ", ", round(psev_q75,  2), "]")
  ) |>
  arrange(dist_l2) |>
  select(id, dist_l2, sumrv_iqr, psev_iqr) |>
  gt() |>
  cols_label(
    id        = "ID",
    dist_l2   = "L2 Distance",
    sumrv_iqr = "sumrv — Median [IQR]",
    psev_iqr  = "psev — Median [IQR]"
  ) |>
  tab_header(title = "Posterior summaries: sumrv and psev by individual")

# PDF — requires webshot2 + chromium
gtsave(tbl_gt, "posterior_summary.pdf")

# RTF
gtsave(tbl_gt, "posterior_summary.rtf")

print(tbl_gt)
