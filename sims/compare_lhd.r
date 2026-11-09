library(ggplot2)
library(dplyr)
library(patchwork)

# Helper Functions

alg_density_plots <- function(sim_results, alpha_value = NULL) {
  sim_config <- sim_results |>
    distinct(n, m)

  plot_list <- list()

  for (i in seq_len(nrow(sim_config))) {
    config_dat <- sim_config[i, ]

    if(!is.null(alpha_value)){
      dat <- sim_results |>
      filter(n == config_dat$n, m == config_dat$m, alpha == alpha_value)
    }
    else{
      dat <- sim_results |> 
        filter(n == config_dat$n, m == config_dat$m)
    }

    p <- ggplot(dat, aes(x = measure, fill = n_chains, group = n_chains)) +
      geom_density(alpha = 0.3) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.001))+
      labs(
        title = sprintf("(%s, %s)", config_dat$n, config_dat$m),
        x = "Measure",
        y = "Density"
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    plot_list[[i]] <- p
    }
  wrap_plots(plot_list) + plot_layout(guides = 'collect')
}

compare_methods_lhd <- function(sim_results, colors = NULL, digits = 0.001) {
  sim_config <- sim_results |>
    distinct(n, m)

  plot_list <- list()

  for (i in seq_len(nrow(sim_config))) {
    config_dat <- sim_config[i, ]

    dat <- sim_results |>
      filter(n == config_dat$n, m == config_dat$m)

    if (!is.null(colors)) {
      lvls <- unique(as.character(dat$alg))
      vals <- colors[lvls]
      fill_scale <- scale_fill_manual(
        values = vals,
        breaks = lvls,
        drop = FALSE
      )
    } else {
      fill_scale <- scale_fill_discrete()
    }

    p <- ggplot(dat, aes(x = measure, fill = alg, group = alg)) +
      geom_density(alpha = 0.4) +
      scale_x_continuous(labels = scales::label_number(accuracy = digits, drop0trailing = TRUE))+
      labs(
        title = sprintf("(%s, %s)", config_dat$n, config_dat$m),
        x = "Measure",
        y = "Density",
        fill = ""
      ) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      fill_scale

    plot_list[[i]] <- p
  }
  n_cols <- dplyr::n_distinct(sim_config$m)
  wrap_plots(plot_list, ncol = n_cols, byrow = FALSE) + plot_layout(guides = 'collect') &
    theme(legend.position = "bottom")
}

# compare PT-maximinLHD (alpha, chains)
pt_maximin_dat <- readRDS(here::here("sims","data","PT_maximinLHD_sim.rds")) |> 
  mutate(across(c(alpha,n_chains),as.factor))

pt_maximin_dat |>
  group_by(n, m, alpha, n_chains) |>
  summarise(
    mean_measure = mean(measure),
    sd_measure = sd(measure),
    lower_bound = mean_measure - sd_measure
  ) |>
  ungroup() |>
  group_by(n, m) |>
  arrange(lower_bound, .by_group = TRUE) |> 
  slice_head(n=2)


alg_density_plots(pt_maximin_dat,alpha_value = 0.25)
alg_density_plots(pt_maximin_dat,alpha_value = 0.5)
alg_density_plots(pt_maximin_dat,alpha_value = 0.7)
alg_density_plots(pt_maximin_dat,alpha_value = 0.9)


# compare PT-maxproLHD (alpha, chains)
pt_maxpro_dat <- readRDS(here::here("sims","data","PT_maxproLHD_sim.rds")) |> 
  mutate(across(c(alpha,n_chains),as.factor))

pt_maxpro_dat |>
  group_by(n, m, alpha, n_chains) |>
  summarise(
    mean_measure = mean(measure),
    sd_measure = sd(measure),
    lower_bound = mean_measure - sd_measure
  ) |>
  ungroup() |>
  group_by(n, m) |>
  arrange(lower_bound, .by_group = TRUE) |> 
  slice_head(n=2)

alg_density_plots(pt_maxpro_dat,alpha_value = 0.25)
alg_density_plots(pt_maxpro_dat,alpha_value = 0.5)
alg_density_plots(pt_maxpro_dat,alpha_value = 0.7)
alg_density_plots(pt_maxpro_dat,alpha_value = 0.9)

# Compare PT with MaxPro R package

pt_maxpro_default <- pt_maxpro_dat |> 
  filter(alpha == 0.9, n_chains == 7) |> 
  select(n,m,measure) |> 
  mutate(alg = "PT")

maxpro_default <- readRDS(here::here("sims","data","maxproLHD_sim.rds")) |> 
  select(n,m,measure) |> 
  mutate(alg = "MaxPro")

maxpro_compare <- rbind(pt_maxpro_default, maxpro_default)

p1 <- compare_methods_lhd(maxpro_compare, 
  colors = c(PT = "#F8766D", MaxPro = "#00BFC4"),
  digits = 0.1) &
  labs(x = expression(Psi["measure"]))

ggsave(
  filename = here::here("sims", "figs", "pt_maxpro_compare.png"),
  plot = p1, height = 7,width = 9, units = "in")


# Compare PT with SLHD R package

pt_maximin_default <- pt_maximin_dat |> 
  filter(alpha == 0.25, n_chains == 7) |> 
  select(n,m,measure) |> 
  mutate(alg = "PT")

slhd_default <- readRDS(here::here("sims","data","SLHD_maximin_sim.rds")) |> 
  select(n,m,measure) |> 
  mutate(alg = "SLHD")

slhd_compare <- rbind(pt_maximin_default, slhd_default)

p2 <- compare_methods_lhd(slhd_compare, 
  colors = c(PT = "#F8766D", SLHD = "#00BFC4"),
  digits = 0.00001) &
  labs(x = expression(phi[p] ~ "measure"))
p2

ggsave(
  filename = here::here("sims", "figs", "pt_slhd_compare.png"),
  plot = p2, height = 7,width = 9, units = "in")
