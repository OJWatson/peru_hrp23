
run_sim_setup <- function(l) {

  plaf <- matrix(
    c(rep(c(0.323, 0.664), l$years)),
    ncol=2,
    byrow=TRUE
  )

  dl <- magenta::drug_list_create(resistance_flag = FALSE,
                                  number_of_resistance_loci = 2,
                                  artemisinin_loci = c(0),
                                  cost_of_resistance = c(1,1),
                                  absolute_fitness_cost_flag = TRUE,
                                  epistatic_logic = NULL,
                                  number_of_drugs = 1,
                                  drugs = list(magenta:::drug_create_default_no_resistance()),
                                  mft_flag = FALSE,
                                  temporal_cycling = -1,
                                  sequential_cycling = -1,
                                  sequential_update = 3,
                                  drug_choice = 0,
                                  partner_drug_ratios = 1)


  magenta::pipeline(EIR = l$EIR,
           years = l$years,
           N = l$N,
           spatial_type = NULL,
           itn_cov = 0,
           num_loci = 2,
           update_save = FALSE,
           update_length = 30,
           human_update_save = FALSE,
           summary_saves_only = TRUE,
           genetics_df_without_summarising = TRUE,
           save_lineages = FALSE,
           full_save = TRUE,
           spatial_incidence_matrix = c(rep(1,l$years)),
           spatial_mosquitoFOI_matrix = c(rep(1,l$years)),
           plaf = plaf,
           drug_list = dl,
           seed = l$seed)

}

run_sim_continue <- function(l) {

  plaf <- matrix(
      rep(c(0,0), l$years),
    ncol=2,
    byrow=TRUE
  )

  dl <- magenta::drug_list_create(resistance_flag = FALSE,
                                  number_of_resistance_loci = 2,
                                  artemisinin_loci = c(0),
                                  cost_of_resistance = c(l$fitness,l$fitness),
                                  absolute_fitness_cost_flag = TRUE,
                                  epistatic_logic = NULL,
                                  number_of_drugs = 1,
                                  drugs = list(magenta:::drug_create_default_no_resistance()),
                                  mft_flag = FALSE,
                                  temporal_cycling = -1,
                                  sequential_cycling = -1,
                                  sequential_update = 3,
                                  drug_choice = 0,
                                  partner_drug_ratios = 1)


  magenta::pipeline(EIR = l$EIR,
                    years = l$years,
                    N = l$N,
                    spatial_type = NULL,
                    itn_cov = c(rep(l$itn_cov, l$itn_years),
                                rep(l$itn_cov2, l$itn_years_2)),
                    num_loci = 2,
                    update_save = TRUE,
                    update_length = 30,
                    human_update_save = TRUE,
                    summary_saves_only = TRUE,
                    genetics_df_without_summarising = TRUE,
                    save_lineages = TRUE,
                    full_save = FALSE,
                    spatial_incidence_matrix = c(rep(0,l$years)),
                    spatial_mosquitoFOI_matrix = c(rep(0,l$years)),
                    saved_state_path = l$saved_state_path,
                    plaf = plaf,
                    drug_list = dl,
                    seed = l$seed)

}

