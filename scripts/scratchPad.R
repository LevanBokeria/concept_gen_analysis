







############## EFFECT SIZES within vs between for Rob #########################
within_stats <- results_table_qc_pass_ptp_analyzed %>%
        group_by(arrangement_pairs, congruency) %>%
        get_summary_stats(phase_2_min_phase_1_ses_1_2_perf, type = 'mean_sd')


within_stats$effect_size <- within_stats$mean/within_stats$sd

arr_1_2_between_effect <- (within_stats$mean[2] - within_stats$mean[1]) / 
        ((within_stats$sd[2] + within_stats$sd[1])/2)

arr_3_4_between_effect <- (within_stats$mean[4] - within_stats$mean[3]) / 
        ((within_stats$sd[4] + within_stats$sd[3])/2)
        


# Playing around with BAYESIAN ANOVA tests
# library(tidyverse)
# library(ggpubr)
# library(rstatix)
# library(pacman)
# library(BayesFactor)
# library(readxl)
# 
# 
# # Clean the env
# rm(list=ls())
# 
# # Get the example dataset
# dat=read.table(url('http://pcl.missouri.edu/exp/aovExample.txt'), head=TRUE)
# 
# # Load the results_table_qc_pass_ptp_analyzed
# results_table_qc_pass_ptp_analyzed <-
#         read_excel(
#                 'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
#         ) %>%
#         as_tibble()
# 
# 
# # Summary statistics
# results_table_qc_pass_ptp_analyzed %>%
#         # filter(arr_phase_1_1 == 3 | arr_phase_1_1 == 15) %>%
#         reorder_levels(arr_phase_1_1, order = c('1','14','3','15')) %>%
#         group_by(arr_phase_1_1,congruency,concept_phase_1) %>%
#         get_summary_stats(phase_2_min_phase_1_ses_1_2_perf, type = 'mean_sd')
# 
# 
# # Turn it all into factors
# results_table_qc_pass_ptp_analyzed <- results_table_qc_pass_ptp_analyzed %>%
#         mutate(ptp = as.factor(ptp),
#                congruency = as.factor(congruency),
#                arr_phase_1_1 = as.factor(arr_phase_1_1),
#                concept_phase_1 = as.factor(concept_phase_1))
# 
# # Filter for arr1-2
# arr_1_2_data <- results_table_qc_pass_ptp_analyzed %>% filter(arr_phase_1_1 == 1 | arr_phase_1_1 == 14)
# arr_3_4_data <- results_table_qc_pass_ptp_analyzed %>% filter(arr_phase_1_1 == 3 | arr_phase_1_1 == 15)
# 
# 
# 
# # Do a Bayesian 3-way ANOVA
# 
# 
# bf_arr_1_2 <- anovaBF(phase_2_min_phase_1_ses_1_2_perf~
#                               congruency*
#                               concept_phase_1*
#                               arr_phase_1_1, 
#              data = arr_1_2_data,
#              whichModels="withmain",
#              whichRandom = "ptp", iterations = 100000)

####################### ANOVA PLAY ########################################
# data("PlantGrowth")
# set.seed(1234)
# PlantGrowth %>% sample_n_by(group, size = 1)
# 
# levels(PlantGrowth$group)
# 
# # Summary stats
# PlantGrowth %>%
#         group_by(group) %>%
#         get_summary_stats(weight, type = "mean_sd")
# 
# ggboxplot(PlantGrowth, x = "group", y = "weight")
# 
# PlantGrowth %>% 
#         group_by(group) %>%
#         identify_outliers(weight)
# 
# # Build the linear model
# model  <- lm(weight ~ group, data = PlantGrowth)
# # Create a QQ plot of residuals
# ggqqplot(residuals(model))
# 
# # Compute Shapiro-Wilk test of normality
# shapiro_test(residuals(model))
# 
# # Check normality assumption by groups
# PlantGrowth %>%
#         group_by(group) %>%
#         shapiro_test(weight)
# 
# # Homogeneity of variance
# plot(model, 1)
# 
# 

