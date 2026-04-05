library(parallel)
library(dplyr)
library(ggplot2)
library(patchwork)
library(PTDesign)

# Helper Functions

plot_m <- function(m_val, data) {
  ggplot(
    filter(data, m == m_val),
    aes(
      x = run_sizes,
      y = measure,
      color = factor(n_chains),
      linetype = factor(n_chains),
      shape = factor(n_chains),
      group = interaction(alpha, n_chains)
    )
  ) +
    geom_line() +
    geom_point() +
    facet_wrap(~alpha, labeller = label_bquote(alpha == .(alpha))) +
    scale_color_manual(
      values = c(
        "#0a0a0a", #black
        "#5ca7ed", # blue
        "#D55E00", # red
        "#45e69b", # green
        "#E69F00", # orange
        "#CC79A7", # pink
        "#56B4E9", # light blue
        "#F0E442" # yellow
      )
    ) +
    scale_linetype_manual(
      values = c("solid", "dashed", "twodash", "dotted", "dotdash")
    ) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7)) + # circle, triangle, square, plus, cross
    labs(
      title = paste0("Maximin (m=", m_val, ")"),
      x = "Run Size (n)",
      y = "Minimum Distance",
      color = "Replicas",
      linetype = "Replicas",
      shape = "Replicas"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

compare_methods_plot <- function(m_val, data) {
  ggplot(filter(data, m == m_val),
         aes(x = run_sizes, y = avg_min_l2_dist,
             color = factor(method),
             linetype = factor(method),
             shape = factor(method))) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = c(
     "#0a0a0a", #black
       "#5ca7ed", # blue
        "#D55E00", # red
        "#45e69b", # green
      "#E69F00", # orange
     "#CC79A7",  # pink
      "#56B4E9",# light blue
      "#F0E442" # yellow
    )) +
    scale_linetype_manual(values = c("solid", "dashed", "twodash", "dotted", "dotdash")) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7)) +  # circle, triangle, square, plus, cross
    labs(
      title = paste0("Maximin (m=",m_val,")"),
      x = "Run Size (n)",
      y = "Minimum Distance",
      color = "method",
      linetype = "method",
      shape = "method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

min_zero <- function(x) min(x[x!=0],na.rm = TRUE)

replicate_sims_oofa <- function(row, n_rep = 100, tolerance = 5000) {
  m <- as.numeric(row["m"])
  run_size <- as.numeric(row["run_sizes"])
  alpha <- as.numeric(row["alpha"])
  n_chains <- as.numeric(row["n_chains"])
  
  vals <- replicate(n_rep, {
    x <- maximin_oofa(n=run_size,k=m,num_chains = n_chains, 
      alpha = alpha, no_improvement = tolerance)
    x$measure
  })
  min_zero(vals)
}

# maximin OofA simuation

m <- 8:10
m_multiples <- seq(2,10,2)

single_grid <- expand.grid(m=m, m_multiples = m_multiples) |> 
  mutate(n = m*m_multiples) |> 
  arrange(m,n)

alpha <- c(0.5, 0.7, 0.9)
n_chains <- c(7,11,15)

result <- do.call(
  rbind,
  lapply(1:nrow(single_grid), function(i) {
    cbind(
      single_grid[i, , drop = FALSE],
      expand.grid(alpha = alpha, n_chains = n_chains)
    )
  })
)

reps <- 100
dat <- do.call(rbind, replicate(reps, result, simplify = FALSE))

n_cores <- detectCores() - 1

measures <- mclapply(
  1:nrow(dat), 
  function(i){
    maximin_oofa(n=dat$n[i], k = dat$m[i], 
      num_chains = dat$n_chains[i], alpha = dat$alpha[i],
      no_improvement = 5000)$measure
  },
  mc.cores = n_cores
)

measures <- unlist(measures)
final_dat <- cbind(dat, measure = measures)

saveRDS(final_dat, file = here::here("sims","data","PT_maximinOofA_sim.rds"))

# Figure: OofA(n,m) with m = 8:10, (n = 2m, 4m,..,10m)

dat <- readRDS(file = here::here("sims","data","PT_maximinOofA_sim.rds"))

fig_dat <- dat |>
  filter(alpha %in% c(0.5, 0.7, 0.9), n_chains %in% c(7, 11, 15)) |>
  rename(run_sizes = n) |> 
  filter(if_all(where(is.numeric), ~ !is.infinite(.))) |> 
  group_by(m, run_sizes, alpha, n_chains)

p1  <- plot_m(8, fig_dat)
p2  <- plot_m(9, fig_dat)
p3 <- plot_m(10, fig_dat)

pt_oofa_sim <- p1 + p2 + p3 + 
  plot_layout(guides = "collect", ncol=1) & theme(legend.position = "bottom")

ggsave(here::here("sims","figs","pt_oofa_sim.png"), 
    plot = pt_oofa_sim, width =7, height = 9, units="in")

# Compare Methods

meta_dat <- readRDS(here::here("sims","data","meta_maximin_m8-10.rds"))

pt_dat <- readRDS(file = here::here("sims","data","PT_maximinOofA_sim.rds"))
fig_dat <- pt_dat |>
  filter(alpha %in% c(0.5, 0.7, 0.9), n_chains %in% c(7, 11, 15)) |>
  as_tibble()

meta_dat <- meta_dat |> 
  filter(method %in% c("DE1","PSO1")) |> 
  mutate(method = factor(sub("1$","",method))) |> 
  rename(run_sizes = n, 
    avg_min_l2_dist = value)

pt_dat <- pt_dat |> 
  filter(alpha == 0.7 & n_chains == 11) |> 
  mutate(method = factor("PT")) |> 
  select(m,run_sizes,method,avg_min_l2_dist)

compare_dat1 <- bind_rows(pt_dat, meta_dat)

fig_compare_oofa <- compare_methods_plot(8, compare_dat1)+
compare_methods_plot(9, compare_dat1)+
compare_methods_plot(10, compare_dat1)+
    plot_layout(guides = "collect", ncol=3) & theme(legend.position = "bottom")

ggsave(here::here("sims","figs","compare_methods_oofa.png"), 
    plot = fig_compare_oofa, width =9, height = 5, units="in")
