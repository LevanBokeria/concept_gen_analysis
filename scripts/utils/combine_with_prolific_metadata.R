combine_with_prolific_metadata = function(dataIn,prolific_meta_data){
        
        # Load libraries
        if (!('tidyverse' %in% (.packages()))){
                library(tidyverse) 
        }
        if (!('readxl' %in% (.packages()))){
                library(readxl)
        }        
        
        
        if(missing(dataIn) | missing(prolific_meta_data)){
                
                # Read the excel files
                dataIn <- read_excel(
                        'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/results_table_all_ptp.xlsx') %>%
                        as_tibble()
                
                # long_form_data_all_ptp <- read_excel(
                #         'C:/Users/levan/GitHub/concept_gen_analysis/results/analysis/long_form_data_all_ptp.xlsx',
                #         sheet = 1) %>%
                #         as_tibble()
                
                prolific_meta_data <- 
                        read_excel(
                                r'(C:\Users\levan\GitHub\concept_gen_analysis\data\prolific_meta_data\united_meta_data.xlsx)',
                                sheet = 1) %>%
                        as_tibble()
        }
        
        # First, create the raw columns in the dataIn 
        dataIn$prolific_status <- NA
        dataIn$prolific_started_datetime <- NA
        dataIn$prolific_completed_date_time <- NA
        dataIn$prolific_time_taken <- NA
        dataIn$prolific_time_taken_min <- NA
        dataIn$prolific_age <- NA
        dataIn$prolific_entered_code <- NA
        dataIn$prolific_sex <- NA
        
        # Go through each participant ID in the prolific meta data
        # Find that person in the dataIn, and add the columns. 
        
        # Find the start date/time for each participant in the results table
        for (i in seq(53,nrow(prolific_meta_data))) {
                
                curr_ptp <- prolific_meta_data$participant_id[i]
                
                # Find the row containing this participant in the results table file
                row_idx <- which(str_detect(dataIn$ptp,curr_ptp))  
                
                if (is_empty(row_idx)){
                        
                        # One participant did not copy over the full ID
                        if (curr_ptp == "5df1117111d87901da0dbc3c"){
                                
                                trimmed_name_to_find <- "5df1117111d87901d"
                                
                                row_idx <- which(str_detect(
                                        dataIn$ptp,
                                        trimmed_name_to_find)
                                )
                        } 
                }
                
                # For that row, copy over the datetime
                dataIn$prolific_status[row_idx] <- 
                        prolific_meta_data$status[i]
                dataIn$prolific_started_datetime[row_idx] <-
                        prolific_meta_data$started_datetime[i]
                dataIn$prolific_completed_date_time[row_idx] <-
                        prolific_meta_data$completed_date_time[i]
                dataIn$prolific_time_taken[row_idx] <- 
                        prolific_meta_data$time_taken[i]
                dataIn$prolific_time_taken_min[row_idx] <- 
                        prolific_meta_data$time_taken_min[i]
                dataIn$prolific_age[row_idx] <- 
                        prolific_meta_data$age[i]
                dataIn$prolific_entered_code[row_idx] <- 
                        prolific_meta_data$entered_code[i]
                dataIn$prolific_sex[row_idx] <- 
                        prolific_meta_data$Sex[i]
                
        }
        
        return(dataIn)
}