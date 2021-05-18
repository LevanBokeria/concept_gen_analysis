# This is a wrapper function to load the datasets, transform columns to factors and add new columns, 
# get only qc_pass participants, and gather some of the columns.

load_transform_gather <- function(){
        
        # Load libraries
        if (!('tidyverse' %in% (.packages()))){
               library(tidyverse) 
        }
        if (!('readxl' %in% (.packages()))){
                library(readxl)
        }    
        if (!('rstatix' %in% (.packages()))){
                library(rstatix)
        }        
        
        
        # Source various scripts
        source('./utils/mutate_results_table_columns.R')
        source('./utils/get_only_qc_pass.R')

        # Read the excel files
        results_table_all_ptp_analyzed <- read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_all_ptp_analyzed.csv') %>%
                as_tibble()
        
        long_form_data_all_ptp_analyzed <- read_csv(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/long_form_data_all_ptp_analyzed.csv') %>%
                as_tibble()
        
        prolific_meta_data <- 
                read_excel(
                r'(C:\Users\levan\GitHub\concept_gen_analysis\data\prolific_meta_data\united_meta_data.xlsx)',
                sheet = 1) %>%
                as_tibble()
        
        # Transform columns to factors and add new columns
        results_table_all_ptp_analyzed <- results_table_all_ptp_analyzed %>%
                mutate_results_table_columns()
        
        long_form_data_all_ptp_analyzed <- long_form_data_all_ptp_analyzed %>%
                mutate_results_table_columns()
        
        
        # Get only QC pass participants
        long_form_data_qc_pass_ptp_analyzed <- 
                get_only_qc_pass(long_form_data_all_ptp_analyzed)
        results_table_qc_pass_ptp_analyzed  <- 
                get_only_qc_pass(results_table_all_ptp_analyzed)
        
        results_table_qc_pass_ptp_analyzed_gathered <- 
                results_table_qc_pass_ptp_analyzed %>%
                select(ptp,congruency,experiment,
                       arr_phase_1_1,arr_phase_2_1,
                       arr_phase_1_name,arr_phase_2_name,
                       concept_phase_1,concept_phase_2,
                       phase_1_ses_1_2_perf,phase_2_ses_1_2_perf,
                       phase_1_lower_left_outlier,phase_2_lower_left_outlier) %>%
                pivot_longer(cols=c(phase_1_ses_1_2_perf,phase_2_ses_1_2_perf),
                            names_to='phase',values_to='ses_1_2_perf')
        
        # Return everything
        returnList <- list(
                'long_form_data_all_ptp_analyzed' = 
                        long_form_data_all_ptp_analyzed,
                'long_form_data_qc_pass_ptp_analyzed' = 
                        long_form_data_qc_pass_ptp_analyzed,
                'results_table_all_ptp_analyzed' = results_table_all_ptp_analyzed,
                'results_table_qc_pass_ptp_analyzed' = results_table_qc_pass_ptp_analyzed,
                'results_table_qc_pass_ptp_analyzed_gathered' = 
                        results_table_qc_pass_ptp_analyzed_gathered)
        return(returnList)
}