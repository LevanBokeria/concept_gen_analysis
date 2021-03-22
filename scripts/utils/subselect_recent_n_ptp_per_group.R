# Description

# This script will take in the results_table_qc_pass_ptp_analyzed file, and the 
# meta-data files from prolific. For each participant, it will copy over the 
# experiment start date and time to the results_table_qc_pass_ptp_analyzed. 
# Then, it will group the table by all the conditions, order by date and time,
# and subselect only N participants per subgroup.

# Load all the libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Clear the environment
rm(list=ls())

# Load the results_table_qc_pass_ptp_analyzed
results_table_qc_pass_ptp_analyzed <-
        read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx'
        ) %>%
        as_tibble()
results_table_qc_pass_ptp_analyzed <- results_table_qc_pass_ptp_analyzed %>%
        mutate(congruency=as.factor(congruency),
               concept_phase_1=as.factor(concept_phase_1),
               concept_phase_2=as.factor(concept_phase_2),
               arr_phase_1_1=as.factor(arr_phase_1_1))

# Read the prolific meta-data files
prolific_meta_1 <-
        read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/data/prolific_meta_data/prolific_export_game_of_image_associations.csv'
        ) %>%
        as_tibble()
prolific_meta_2 <-
        read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/data/prolific_meta_data/prolific_export_game_of_finding_image_associations.csv'
        ) %>%
        as_tibble()

# Combine the two prolific meta_data files
prolific_meta_data <- bind_rows(prolific_meta_1,prolific_meta_2, .id = "id")

###############################################################################
###############################################################################

# Add date and time column to the results table
results_table_qc_pass_ptp_analyzed$started_datetime <- as.POSIXct(NA)
results_table_qc_pass_ptp_analyzed$completed_date_time <- as.POSIXct(NA)


# Find the start date/time for each participant in the results table
for (i in seq(1,nrow(prolific_meta_data))) {
        
        curr_ptp <- prolific_meta_data$participant_id[i]
        
        # Find the row containing this participant in the results table file
        row_idx <- which(str_detect(results_table_qc_pass_ptp_analyzed$ptp,curr_ptp))  
        
        if (is_empty(row_idx)){
                
                # One participant did not copy over the full ID
                if (curr_ptp == "5df1117111d87901da0dbc3c"){
                        
                        trimmed_name_to_find <- "5df1117111d87901d"
                        
                        row_idx <- which(str_detect(
                                results_table_qc_pass_ptp_analyzed$ptp,
                                trimmed_name_to_find)
                                )
                } 
        }
        
        
        # For that row, copy over the datetime
        results_table_qc_pass_ptp_analyzed$started_datetime[row_idx] <-
                prolific_meta_data$started_datetime[i]
        results_table_qc_pass_ptp_analyzed$completed_date_time[row_idx] <-
                prolific_meta_data$completed_date_time[i]
        
        }

# Make sure every participant has a date time in the results table


# Separate the first two arrangement data from the second two
results_arr_pairs_1 <- results_table_qc_pass_ptp_analyzed %>%
        subset(arr_phase_1_1 == "1" | arr_phase_1_1 == "14")
results_arr_pairs_2 <- results_table_qc_pass_ptp_analyzed %>%
        subset(arr_phase_1_1 == "3" | arr_phase_1_1 == "15")

# Group the results table
results_arr_pairs_1_balanced <- results_arr_pairs_1 %>%
        group_by(congruency,concept_phase_1,arr_phase_1_1) %>%
        arrange(completed_date_time) %>%
        slice_head(n = 8) %>%
        ungroup()

results_arr_pairs_2_balanced <- results_arr_pairs_2 %>%
        group_by(congruency,concept_phase_1,arr_phase_1_1) %>%
        arrange(completed_date_time) %>% 
        slice_head(n = 15) %>%
        ungroup()

# Combine them now
results_balanced <- bind_rows(results_arr_pairs_1_balanced,
                              results_arr_pairs_2_balanced)


# Write the balanced results as a csv file
results_balanced %>%
        write_csv("../results/analysis/results_table_qc_pass_ptp_analyzed_balanced.csv")
