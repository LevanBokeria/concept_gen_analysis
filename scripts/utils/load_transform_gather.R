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
        
        # Source various scripts
        source('./utils/mutate_results_table_columns.R')
        source('./utils/get_only_qc_pass.R')
        source('./utils/basic_checks.R')
        
        # Read the excel files
        results_table_qc_pass_ptp_analyzed <- read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_qc_pass_ptp_analyzed.xlsx') %>%
                as_tibble()
        
        long_form_data_all_ptp <- read_excel(
                'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/long_form_data_all_ptp.xlsx',
                sheet = 1) %>%
                as_tibble()
        
        prolific_meta_data <- 
                read_excel(
                r'(C:\Users\levan\GitHub\concept_gen_analysis\data\prolific_meta_data\united_meta_data.xlsx)',
                sheet = 1) %>%
                as_tibble()
        
        # Basic checks
        basic_checks(long_form_data_all_ptp)
        
        # Transform columns to factors and add new columns
        results_table_qc_pass_ptp_analyzed <- results_table_qc_pass_ptp_analyzed %>%
                mutate_results_table_columns()
        
        long_form_data_all_ptp <- long_form_data_all_ptp %>%
                mutate_results_table_columns()
        
        
        # Get only QC pass participants
        long_form_data_qc_pass_ptp <- get_only_qc_pass(long_form_data_all_ptp)
        
        
        # Gather, otherwise its not actually in a long format
        results_table_qc_pass_ptp_analyzed_gathered <- 
                results_table_qc_pass_ptp_analyzed %>%
                select(ptp,congruency,arrangement_pairs,arr_phase_1_1,arr_phase_2_1,
                       phase_1_ses_1_2_perf,phase_2_ses_1_2_perf,
                       phase_1_lower_left_outlier,phase_2_lower_left_outlier) %>%
                gather(key='phase',value='ses_1_2_perf',
                       phase_1_ses_1_2_perf,
                       phase_2_ses_1_2_perf)
        
        
        # Return everything
        returnList <- list('long_form_data_all_ptp' = long_form_data_all_ptp,
                           'long_form_data_qc_pass_ptp' = long_form_data_qc_pass_ptp,
                           'results_table_qc_pass_ptp_analyzed' = 
                                   results_table_qc_pass_ptp_analyzed,
                           'results_table_qc_pass_ptp_analyzed_gathered' =
                                   results_table_qc_pass_ptp_analyzed_gathered)
        return(returnList)
}