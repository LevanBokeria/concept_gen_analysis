# Load libraries
library(tidyverse)
library(readxl)
library(ggpubr)
library(rstatix)
library(emmeans)


# Clear the environment
rm(list = ls())

# Load the results_table_qc_pass_ptp_analyzed
# results_table_qc_pass_ptp_analyzed <-
#         read_excel(
#                 'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
#         ) %>%
#         as_tibble()

# Load the BALANCED results_table_qc_pass_ptp_analyzed
results_table_qc_pass_ptp_analyzed <-
        read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed_balanced.csv'
        ) %>%
        as_tibble()

results_table_qc_pass_ptp_analyzed <- results_table_qc_pass_ptp_analyzed %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

# Load the long_form_data_all_ptp
long_form_data_all_ptp <-
        read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/long_form_data_all_ptp.xlsx'
                , sheet = 1
        ) %>%
        as_tibble()

long_form_data_all_ptp <- long_form_data_all_ptp %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

long_form_data_qc_pass_ptp <- subset(long_form_data_all_ptp, 
                                     data_submitted == 1 & 
                                     global_pass_incl_phase_2_fails == 1 & 
                                     debrief_qc_pass == 1 & 
                                     fb_int_qc_pass == 1 & 
                                     phase_1_rt_qc_pass == 1 & 
                                     phase_2_rt_qc_pass == 1)




# ###################### FILTER THE DATA ###################################
# long_form_data_qc_pass_ptp <- 
#         subset(long_form_data_qc_pass_ptp, 
#                                      arr_phase_1_1 == 1 | arr_phase_1_1 == 14)
# results_table_qc_pass_ptp_analyzed <- 
#         subset(results_table_qc_pass_ptp_analyzed, 
#                                              arr_phase_1_1 == 1 | arr_phase_1_1 == 14)


###################### GATHER ##############################################
# Gather the long form data for results_table_qc_pass_ptp_analyzed
results_table_qc_pass_ptp_analyzed_gathered <- 
        results_table_qc_pass_ptp_analyzed %>%
        select(ptp,congruency,phase_1_ses_1_2_perf,phase_2_ses_1_2_perf) %>%
        gather(key='phase',value='ses_1_2_perf',
               phase_1_ses_1_2_perf,
               phase_2_ses_1_2_perf)


####################### START ANOVA ANALYSIS ###################################

levels(results_table_qc_pass_ptp_analyzed$arr_phase_1_1)

# Summary statistics
results_table_qc_pass_ptp_analyzed %>%
        group_by(congruency,arr_phase_1_1) %>%
        get_summary_stats(phase_2_min_phase_1_ses_1_2_perf, type = 'mean_sd')

# Check for outliers
# results_table_qc_pass_ptp_analyzed %>%
#         group_by(congruency) %>%
#         identify_outliers(phase_2_min_phase_1_ses_1_2_perf) %>%
#         View()

# ggqqplot(results_table_qc_pass_ptp_analyzed, "phase_2_min_phase_1_ses_1_2_perf", facet.by = "congruency")



# Add arr1-2 VS arr3-4 as factors
results_table_qc_pass_ptp_analyzed <- 
        results_table_qc_pass_ptp_analyzed %>%
        mutate(arrangement_pairs = "arr_1_2")

results_table_qc_pass_ptp_analyzed$arrangement_pairs[
        results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "3" | 
                results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "15"
] <- "arr_3_4"

results_table_qc_pass_ptp_analyzed <-
results_table_qc_pass_ptp_analyzed %>%
        mutate(arrangement_pairs = as.factor(arrangement_pairs))

# Add arr 1,2,3,4 properly formatted to be sub factors for arrangement_pairs
results_table_qc_pass_ptp_analyzed$arr_idx <- "arr_idx_1"

results_table_qc_pass_ptp_analyzed$arr_idx[
        results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "1"] <- "arr_idx_1"
results_table_qc_pass_ptp_analyzed$arr_idx[
        results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "14"] <- "arr_idx_2"
results_table_qc_pass_ptp_analyzed$arr_idx[
        results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "3"] <- "arr_idx_2"
results_table_qc_pass_ptp_analyzed$arr_idx[
        results_table_qc_pass_ptp_analyzed$arr_phase_1_1 == "15"] <- "arr_idx_1"

results_table_qc_pass_ptp_analyzed <- 
        results_table_qc_pass_ptp_analyzed %>%
        mutate(arr_idx = as.factor(arr_idx))


# Do a 4-way anova, with congruency concept order, arrangement order, and arrangement pair
res.aov_4 <- results_table_qc_pass_ptp_analyzed %>% anova_test(
        phase_2_min_phase_1_ses_1_2_perf ~ 
                congruency*concept_phase_1*arr_idx*arrangement_pairs)
res.aov_4

res.aov_3_arr_1_2 <- results_table_qc_pass_ptp_analyzed %>% 
        filter(arrangement_pairs == "arr_1_2") %>%
        anova_test(
        phase_2_min_phase_1_ses_1_2_perf ~ 
                congruency*concept_phase_1*arr_phase_1_1)
res.aov_3_arr_1_2

res.aov_3_arr_3_4 <- results_table_qc_pass_ptp_analyzed %>% 
        filter(arrangement_pairs == "arr_3_4") %>%
        anova_test(
                phase_2_min_phase_1_ses_1_2_perf ~ 
                        congruency*concept_phase_1*arr_phase_1_1)
res.aov_3_arr_3_4

bxp <- ggboxplot(
        results_table_qc_pass_ptp_analyzed, x = "arr_phase_1_1", y = "phase_2_min_phase_1_ses_1_2_perf",
        color = "congruency", palette = "jco"
)
bxp
