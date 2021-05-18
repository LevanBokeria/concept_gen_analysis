mutate_results_table_columns = function(dataIn){
        
        # print('starting the script')
        
        # For both long_form and results_table files
        dataIn <- dataIn %>%
                mutate(across(c(ptp,progress_state,data_submitted,
                                concept_phase_1,concept_phase_2,
                                arr_phase_1_1,arr_phase_2_1,
                                congruency,
                                global_pass_incl_phase_2_fails,
                                arr_phase_1_name,arr_phase_2_name,
                                phase_1_lower_left_outlier,
                                phase_2_lower_left_outlier,
                                experiment,
                                basic_data_checks_pass),
                              as.factor))
        
        # If its the results table, then 'status' exists. Mutate the prolific 
        # data variables
        if ('status' %in% colnames(dataIn)){
          dataIn <- dataIn %>%
                  mutate(across(c(status,
                                  entered_code,
                                  Sex),
                                as.factor))
          
          # Reorder levels for various things
          dataIn <-
                  dataIn %>%
                  reorder_levels(arr_phase_1_1, order = c(1,14,3,15,NaN)) %>%
                  reorder_levels(congruency, order = c(0,1,NaN)) %>%
                  reorder_levels(arr_phase_1_name, order = c(1,2,3,4,NaN)) %>%
                  reorder_levels(arr_phase_2_name, order = c(1,2,3,4,NaN))
          
        }
        
        
        # If prompt img name exists, mutate that too
        if ('prompt_img_name' %in% colnames(dataIn)){
                dataIn <- dataIn %>% 
                        mutate(prompt_img_name=as.factor(prompt_img_name))
                
                # Reorder levels for various things
                dataIn <-
                        dataIn %>%
                        reorder_levels(arr_phase_1_1, order = c(1,14,3,15)) %>%
                        reorder_levels(congruency, order = c(0,1)) %>%
                        reorder_levels(arr_phase_1_name, order = c(1,2,3,4)) %>%
                        reorder_levels(arr_phase_2_name, order = c(1,2,3,4))
        }

        return(dataIn)
}