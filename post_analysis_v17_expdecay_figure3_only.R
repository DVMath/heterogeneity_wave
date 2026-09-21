# ---------------------------------------------------------------------------
# Fast Figure 3 (sir_extreme_individuals) re-plot from a saved checkpoint.
# Only valid for pure styling/labeling changes (theme, fonts, titles, legend,
# colour-scale fixes) -- if anything upstream changes (data, model,
# archetype selection, simulation logic), re-run the full
# post_analysis_v17_expdecay_concurrent.R instead; this script does not
# recompute any of that.
# Requires checkpoint_before_fig3.RData, saved by that script right after its
# two long simulation loops finish (see the "Checkpoint" comment there).
# ---------------------------------------------------------------------------
library(tidyverse)
library(rstan)
library(patchwork)
library(cowplot)
library(glue)

stopifnot(file.exists("checkpoint_before_fig3.RData"))
load("checkpoint_before_fig3.RData")
cat("Loaded checkpoint_before_fig3.RData\n")

cat("Saved checkpoint_before_fig3.RData\n")

my_theme_pub_y <- theme_classic(base_size = 20, base_family = "serif") +
  theme(
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 20),
    axis.title       = element_text(size = 22),
    plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(colour = "grey94", linewidth = 0.3),
    strip.text       = element_text(size = 19, face = "bold"),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 20),
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
  scale_colour_manual(values = scenario_colours_5,
                      labels = scenario_labels_5, name = NULL) +
  scale_fill_manual(  values = scenario_colours_5, guide = "none") +
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
  scale_colour_manual(values = scenario_colours_5,
                      labels = scenario_labels_5, name = NULL) +
  scale_fill_manual(  values = scenario_colours_5, guide = "none") +
  scale_y_continuous(labels = scales::comma, limits = deaths_lim) +
  facet_wrap(~ scenario, nrow = 1, labeller = sc_labeller_5) +
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
       title = "Population level") +
  my_theme_pub_y +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())

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
          legend.text      = element_text(size = 20),
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
      plot.tag          = element_text(size = 16, face = "bold", family = "serif"),
      plot.tag.position = "topleft"
    )
  )

p_final

ggsave(outpath("sir_extreme_individuals.pdf"),  plot = p_final,
       width = 20, height = 12, device = cairo_pdf)
ggsave(outpath("sir_extreme_individuals.tiff"), plot = p_final,
       width = 20, height = 12, dpi = 600, compression = "lzw")

cat("Saved sir_extreme_individuals.pdf and .tiff\n")
