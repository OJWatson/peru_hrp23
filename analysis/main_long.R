# Start with installing needed packages
install.packages(
  "hipercow",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org")
)

install.packages(
  "rrq",
  repos = c("https://mrc-ide.r-universe.dev", "https://cloud.r-project.org")
)


# Now load up a new driver
library(hipercow)
library(tidyverse)
hipercow_init()

#And set up configuration
hipercow_configure(driver = "dide-windows")
hipercow_configuration()
windows_authenticate()

# Check basic works
id <- task_create_expr(sessionInfo())
task_status(id)

# setup packages
hipercow_provision()

# check magenta runs
id <- task_create_expr(magenta::pipeline(EIR = 10,
                                         years = 1, N = 1000,
                                         spatial_type = NULL,
                                         num_loci = 2,
                                         update_save = TRUE,
                                         update_length = 30,
                                         human_update_save = TRUE,
                                         summary_saves_only = TRUE,
                                         genetics_df_without_summarising = TRUE,
                                         save_lineages = TRUE,
                                         full_save = FALSE,
                                         seed = 123456789L))
task_status(id)
task_log_show(id)

# setup environment
source("analysis/funcs.R")
hipercow_environment_create(sources = c("analysis/funcs.R"), overwrite = TRUE)

# HELPER AND PLOTTING FUNCTIONS --------------

# helper function for getting task paths
get_result_path <- function(id, follow = TRUE) {
  root <- hipercow:::hipercow_root(NULL)
  id <- hipercow:::check_task_id(id, "task_result", TRUE, call = rlang::current_env())
  if (follow) {
    id <- hipercow:::follow_retry_map(id, root)
  }
  path <- hipercow:::path_task(root$path$tasks, id)
  path_result <- file.path(path, hipercow:::RESULT)
  path_result
}

# get cases
get_cases <- function(x, N = 550000) {

  return(
    data.frame(
      cases_per_cap = unlist(lapply(seq_len(length(x)-1), function(i) {sum(unlist(x[[i]][3:5]))/N*1000})),
      pcr_prev = unlist(lapply(seq_len(length(x)-1), function(i) {sum(unlist(x[[i]]$pcr_prev))})),
      l0 = unlist(lapply(seq_len(length(x)-1), function(i) {x[[i]]$lineage[1]})),
      l1 = unlist(lapply(seq_len(length(x)-1), function(i) {x[[i]]$lineage[2]})),
      l2 = unlist(lapply(seq_len(length(x)-1), function(i) {x[[i]]$lineage[3]})),
      l3 = unlist(lapply(seq_len(length(x)-1), function(i) {x[[i]]$lineage[4]})),
      t = seq(0, (length(x)-2)*30, 30)
    )
  )
}

# create case df
create_case_df <- function(res, parl) {
  do.call(rbind, lapply(seq_along(res), function(x) {
    df <- get_cases(res[[x]], parl[[x]]$N[1])
    df$N <- parl[[x]]$N[1]
    df$seed <- parl[[x]]$seed[1]
    return(df)
  }))
}

# create dfp
create_dfp <- function(case_df) {
  case_df %>%
    filter(pcr_prev  > 0) %>%
    mutate(t = t/365) %>% filter(t > 0) %>%
    mutate(n = l0+l1+l2+l3) %>%
    mutate(WT = l0/n, hrp2d = l1/n, hrp3d = l2/n, hrp23d = l3/n) %>%
    pivot_longer(WT:hrp23d) %>% mutate(Genotype = name)
}

# gen plot
gen_plot <- function(dfp) {
  gg <- dfp %>%
    mutate(t = 2006+t) %>%
    mutate(N = factor(paste("N =", N), levels =  paste("N =", seq(25000, 200000, 25000)))) %>%
    ggplot(aes(t, value, color = Genotype, group = interaction(seed, Genotype))) +
    geom_line(alpha = 0.1) +
    geom_line(aes(t, value, color = Genotype), data = dfp %>% group_by(t, Genotype) %>% summarise(value = mean(value,na.rm=TRUE)) %>%
                mutate(t = 2006+t), inherit.aes = FALSE, linewidth = 2) +
    theme_bw(base_size = 10) +
    theme(axis.line = element_line(), axis.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_vline(xintercept = c(2006, 2010)) +
    scale_x_continuous(breaks = c(2006,2010,2014,2018)) +
    facet_grid(N~Genotype) +
    ylab("\nGenotype Frequency") +
    xlab("Year") +
    MetBrewer::scale_color_met_d("Egypt")
  print(gg)
  invisible(gg)
}

# cases plot
cases_plot <- function(dfp){
  gg <- dfp %>%
    mutate(t = 2006+t) %>%
    mutate(N = factor(paste("N =", N), levels =  paste("N =", seq(25000, 200000, 25000)))) %>%
    ggplot(aes(t, cases_per_cap, group = seed)) +
    geom_line(alpha = 0.1) +
    geom_vline(xintercept = c(2006, 2010)) +
    geom_line(aes(t, cases_per_cap), data = dfp %>% group_by(t, Genotype) %>% summarise(cases_per_cap = median(cases_per_cap)) %>%
                mutate(t = 2006+t), inherit.aes = FALSE, linewidth = 2) +
    theme_bw(base_size = 10) +
    theme(axis.line = element_line(), axis.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Monthly Cases / 1000") +
    xlab("Year")+
    scale_x_continuous(breaks = c(2006,2010,2014,2018)) +
    scale_y_sqrt(breaks = c(0,0.25, 1,2)) +
    facet_grid(N~.)
  print(gg)
  invisible(gg)
}

# FINAL # -----------------

## final setup ------------------
pardf6 <- data.frame("EIR" = 0.101,
                     "years" = 30,
                     "N" = seq(25000, 200000, 25000),
                     "rep" = sort(rep(1:500, 8)))
pardf6$seed <- seq_len(nrow(pardf6))
parl6 <- split(pardf6, pardf6$seed)

# bulk submit setup
# N_setup_fin_test_full
bundle6 <- task_create_bulk_call(run_sim_setup, data = parl6, bundle_name = "N_setup_fin_long_full")
while(any(hipercow::task_status(bundle6$ids) != "success")) {
  Sys.sleep(60)
}
 res6 <- lapply(bundle6$ids, hipercow::task_result)

# get the results paths and modify
for(i in seq_along(bundle6$ids)) {
  message(i)
  out <- hipercow::task_result(bundle6$ids[i])
  if(all(out$population_List$Infection_States %in% c(0,5)) ||
     all(lengths(out$scourge_List$Mosquito_Oocyst_barcode_male_vectors)==0) ||
     all(lengths(out$scourge_List$Mosquito_Oocyst_barcode_female_vectors)==0) ||
     all(lengths(out$populations_event_and_strains_List$Strain_barcode_vectors)==0) ||
     all(all(lengths(out$populations_event_and_strains_List$Infection_barcode_realisation_vectors)==0) )) {
    out <- hipercow::task_result(bundle6$ids[which(pardf6$N == pardf6$N[i])[1]])
  }
  out <- magenta::update_saved_state_barcode_plaf(out, c(0.323, 0.664))
  saveRDS(out, file.path("analysis/results/setup/", paste0(bundle6$ids[i], ".rds")))
  gc()
}
paths <- file.path("analysis/results/setup/", paste0(bundle6$ids, ".rds"))


pardf7 <- data.frame("EIR" = 0.101,
                     "years" = 4 + 8,
                     "itn_cov" = 0.2,
                     "itn_cov2" = 0,
                     "itn_years" = 4,
                     "itn_years_2" = 8,
                     "N" = seq(25000, 200000, 25000),
                     "rep" = sort(rep(1:100, 8)),
                     "saved_state_path" = paths)
pardf7$seed <- seq_len(nrow(pardf7))
parl7 <- split(pardf7, pardf7$seed)

# bulk submit continue
bundle7 <- task_create_bulk_call(run_sim_continue, data = parl7, bundle_name = "N_continue_fin_long_full")
hipercow::task_status(bundle7$ids) %>% table
while(any(hipercow::task_status(bundle7$ids) != "success")) {
  Sys.sleep(60)
}
res7 <- lapply(bundle7$ids, hipercow::task_result)
saveRDS(res7, "analysis/results/sim_N_long.rds")

# Get results
case_df7 <- create_case_df(res7, parl7)
dfp7 <- create_dfp(case_df7)

# Do without itn changes
pardf8 <- data.frame("EIR" = 0.101,
                     "years" = 4 + 8,
                     "itn_cov" = 0.0,
                     "itn_cov2" = 0.0,
                     "itn_years" = 4,
                     "itn_years_2" = 8,
                     "N" = seq(25000, 200000, 25000),
                     "rep" = sort(rep(1:100, 8)),
                     "saved_state_path" = paths)
pardf8$seed <- seq_len(nrow(pardf8))
parl8 <- split(pardf8, pardf8$seed)

bundle8 <- task_create_bulk_call(run_sim_continue, data = parl8, bundle_name = "N_continue_no_itn_long_test_full")
hipercow::task_status(bundle8$ids) %>% table
while(any(hipercow::task_status(bundle8$ids) != "success")) {
  Sys.sleep(60)
}
res8 <- lapply(bundle8$ids, hipercow::task_result)
saveRDS(res8, "analysis/results/sim_N_no_itn_long.rds")

# Get results
case_df8 <- create_case_df(res8, parl8)
dfp8 <- create_dfp(case_df8)

## final plotting ------------------

cases_gg_7 <- cases_plot(dfp7)
gen_gg_7 <- gen_plot(dfp7)

cases_gg_8 <- cases_plot(dfp8)
gen_gg_8 <- gen_plot(dfp8)

# now probability of event

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

dfp_all_7 <- dfp7 %>% mutate(t = 2006+t) %>% filter(t >= 2016) %>%
  group_by(seed, N) %>%
  summarise(m = mean(value[name == "hrp23d"],na.rm=TRUE)>0.75) %>% group_by(N) %>% summarise(p = sum(m)/n(), n = n()) %>%
  mutate(low = Hmisc::binconf(.data$p*50, 50)[,2]) %>%
  mutate(high = Hmisc::binconf(.data$p*50, 50)[,3]) %>%
  mutate(pamafro = TRUE)

dfp_all_8 <- dfp8 %>% mutate(t = 2006+t) %>% filter(t >= 2016) %>%
  group_by(seed, N) %>%
  summarise(m = mean(value[name == "hrp23d"],na.rm=TRUE)>0.75) %>% group_by(N) %>% summarise(p = sum(m)/n(), n = n()) %>%
  mutate(low = Hmisc::binconf(.data$p*50, 50)[,2]) %>%
  mutate(high = Hmisc::binconf(.data$p*50, 50)[,3]) %>%
  mutate(pamafro = FALSE)

dfp_ind_7 <- dfp7 %>% mutate(t = 2006+t) %>% filter(t >= 2016) %>%
  filter(name == "hrp23d") %>%
  group_by(seed, N) %>%
  summarise(p = mean(value[name == "hrp23d"],na.rm=TRUE)>0.75) %>%
  mutate(pamafro = TRUE)

dfp_ind_8 <- dfp8 %>% mutate(t = 2006+t) %>% filter(t >= 2016) %>%
  filter(name == "hrp23d") %>%
  group_by(seed, N) %>%
  summarise(p = mean(value[name == "hrp23d"],na.rm=TRUE)>0.75) %>%
  mutate(pamafro = FALSE)

prob_gg <- rbind(dfp_all_7, dfp_all_8) %>%
  ggplot(aes((N), p, ymin=low, ymax=high, color = pamafro)) +
  geom_errorbar(position = position_dodge(width = 5000)) +
  geom_point(position = position_dodge(width = 5000), size = 2) +
  binomial_smooth(data = rbind(dfp_ind_7, dfp_ind_8), mapping = aes((N), as.integer(p), color = pamafro, fill = pamafro),
                  alpha = 0.2, inherit.aes = FALSE) +
  ggpubr::theme_pubr(base_size = 10) +
  scale_x_continuous(breaks = seq(25000,200000,25000)) +
  theme(axis.line = element_line(), axis.title = element_text(size = 12), panel.grid.major = element_line(color = "grey"),
        legend.position = "right") +
  xlab("\nHuman Population Size") +
  MetBrewer::scale_fill_met_d(palette_name = "Austria", direction = -1, name = "PAMAFRO") +
  MetBrewer::scale_color_met_d(palette_name = "Austria", direction = -1, name = "PAMAFRO") +
  ylab("Probability of hrp2-/hrp3- frequency >75% in 2016-2018\n")
prob_gg

save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      root = "analysis/plots",
                      dpi = 600) {

  dir.create(root, showWarnings = FALSE)
  fig_path <- function(name) {paste0(root, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     base_width = width, dpi = dpi)

  pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height, )
  print(fig)
  dev.off()


}

save_figs("supp_traces_long", cowplot::plot_grid(cases_gg_7, gen_gg_7, rel_widths = c(1,2.5)), width = 6.3*1.55, height = 8*1.25, dpi = 300)
save_figs("supp_traces_no_pamafro_long", cowplot::plot_grid(cases_gg_8, gen_gg_8, rel_widths = c(1,2.5)), width = 6.3*1.55, height = 8*1.25, dpi = 300)
save_figs("main_prob_long", prob_gg, width = 6.3*1.5, height = 3.5*1.5)

# Reproducibility script ---
dir.create("analysis/env", showWarnings = FALSE)
write.csv(installed.packages("hipercow/lib/windows/4.5.1/"),
          "analysis/env/installed_pacakages.csv")
