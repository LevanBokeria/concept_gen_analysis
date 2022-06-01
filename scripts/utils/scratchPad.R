
rm(list=ls())

source('./utils/load_all_libraries.R')


# How many participants do I still need to run? ##############################

source('./utils/load_transform_gather.R')

outList <- load_transform_gather()

long_form_data_all_ptp_analyzed <- 
        outList$long_form_data_all_ptp_analyzed
long_form_data_qc_pass_ptp_analyzed <- 
        outList$long_form_data_qc_pass_ptp_analyzed

results_table_all_ptp_analyzed <- 
        outList$results_table_all_ptp_analyzed
results_table_qc_pass_ptp_analyzed <- 
        outList$results_table_qc_pass_ptp_analyzed
results_table_qc_pass_ptp_analyzed_gathered <- 
        outList$results_table_qc_pass_ptp_analyzed_gathered

rm(outList)

## Count the number of participants per subgroup --------------------------


results_table_qc_pass_ptp_analyzed %>%
        filter(experiment == 1) %>%
        group_by(congruency,
                 arr_phase_1_name,
                 concept_phase_1) %>%
        summarise(n = n())

results_table_qc_pass_ptp_analyzed %>%
        filter(experiment == 1) %>%
        group_by(congruency,
                 arr_phase_1_name) %>%
        summarise(n = n())

results_table_qc_pass_ptp_analyzed %>%
        filter(experiment == 1) %>%
        group_by(congruency,
                 concept_phase_1) %>%
        summarise(n = n())

results_table_qc_pass_ptp_analyzed %>%
        filter(experiment == 1) %>%
        group_by(concept_phase_1,
                 arr_phase_1_name) %>%
        summarise(n = n())
