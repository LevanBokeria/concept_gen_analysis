mutate_results_table_columns = function(dataIn){
        
        print('starting the script')
        
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
                                status,
                                age,
                                entered_code,
                                Sex,
                                basic_data_checks_pass),
                              as.factor))
        
        # If prompt img name exists, mutate that too
        if ('prompt_img_name' %in% colnames(dataIn)){
                dataIn <- dataIn %>% 
                        mutate(prompt_img_name=as.factor(prompt_img_name))
        }
        
        # # Add a factor describing whether the arrangement has an outlier in the lower left of the space
        # dataIn$phase_1_lower_left_outlier <- TRUE
        # 
        # dataIn$phase_1_lower_left_outlier[
        #         dataIn$arr_phase_1_1 == "1"] <- TRUE
        # dataIn$phase_1_lower_left_outlier[
        #         dataIn$arr_phase_1_1 == "14"] <- FALSE
        # dataIn$phase_1_lower_left_outlier[
        #         dataIn$arr_phase_1_1 == "3"] <- FALSE
        # dataIn$phase_1_lower_left_outlier[
        #         dataIn$arr_phase_1_1 == "15"] <- TRUE
        # 
        # # Same for phase 2
        # dataIn$phase_2_lower_left_outlier <- TRUE
        # 
        # dataIn$phase_2_lower_left_outlier[
        #         dataIn$arr_phase_2_1 == "1"] <- TRUE
        # dataIn$phase_2_lower_left_outlier[
        #         dataIn$arr_phase_2_1 == "14"] <- FALSE
        # dataIn$phase_2_lower_left_outlier[
        #         dataIn$arr_phase_2_1 == "3"] <- FALSE
        # dataIn$phase_2_lower_left_outlier[
        #         dataIn$arr_phase_2_1 == "15"] <- TRUE        
        # 
        
        # Reorder levels for various things
        dataIn <-
                dataIn %>%
                reorder_levels(arr_phase_1_1, order = c(1,14,3,15)) %>%
                reorder_levels(congruency, order = c(0,1))
        
        return(dataIn)
}