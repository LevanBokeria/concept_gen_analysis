mutate_results_table_columns = function(dataIn){
        
        dataIn <- dataIn %>% 
                mutate(congruency=as.factor(as.numeric(congruency)),
                       concept_phase_1=as.factor(concept_phase_1),
                       concept_phase_2=as.factor(concept_phase_2),
                       arr_phase_1_1=as.factor(arr_phase_1_1),
                       arr_phase_2_1=as.factor(arr_phase_2_1))
        
        
        # If ptp exists mutate that too
        if ('ptp' %in% colnames(dataIn)){
                dataIn <- dataIn %>%
                        mutate(ptp=as.factor(ptp))
        }
        
        # If prompt img name exists, mutate that too
        if ('prompt_img_name' %in% colnames(dataIn)){
                dataIn <- dataIn %>% 
                        mutate(prompt_img_name=as.factor(prompt_img_name))
        }        
        
        
        
        # Add arr1-2 VS arr3-4 as factors
        dataIn <- 
                dataIn %>%
                mutate(arrangement_pairs = "arr_1_2")
        
        dataIn$arrangement_pairs[
                dataIn$arr_phase_1_1 == "3" | 
                        dataIn$arr_phase_1_1 == "15"
        ] <- "arr_3_4"
        
        dataIn <-
                dataIn %>%
                mutate(arrangement_pairs = as.factor(arrangement_pairs))
        
        # Add a factor describing whether the arrangement has an outlier in the lower left of the space
        dataIn$phase_1_lower_left_outlier <- TRUE
        
        dataIn$phase_1_lower_left_outlier[
                dataIn$arr_phase_1_1 == "1"] <- TRUE
        dataIn$phase_1_lower_left_outlier[
                dataIn$arr_phase_1_1 == "14"] <- FALSE
        dataIn$phase_1_lower_left_outlier[
                dataIn$arr_phase_1_1 == "3"] <- FALSE
        dataIn$phase_1_lower_left_outlier[
                dataIn$arr_phase_1_1 == "15"] <- TRUE
        
        # Same for phase 2
        dataIn$phase_2_lower_left_outlier <- TRUE
        
        dataIn$phase_2_lower_left_outlier[
                dataIn$arr_phase_2_1 == "1"] <- TRUE
        dataIn$phase_2_lower_left_outlier[
                dataIn$arr_phase_2_1 == "14"] <- FALSE
        dataIn$phase_2_lower_left_outlier[
                dataIn$arr_phase_2_1 == "3"] <- FALSE
        dataIn$phase_2_lower_left_outlier[
                dataIn$arr_phase_2_1 == "15"] <- TRUE        
        
        dataIn <- 
                dataIn %>% 
                mutate(phase_1_lower_left_outlier = as.factor(phase_1_lower_left_outlier),
                       phase_2_lower_left_outlier = as.factor(phase_2_lower_left_outlier))
        
        # Reorder levels for various things
        dataIn <-
                dataIn %>%
                reorder_levels(arr_phase_1_1, order = c(1,14,3,15)) %>%
                reorder_levels(congruency, order = c(0,1))
        
        return(dataIn)
}