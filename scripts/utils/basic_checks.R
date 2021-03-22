# Description:

# Performs some basic checks on the dataset provided.

basic_checks <- function(dataIn){
        
        # 1. Check that the column "global_pass_including_phase_2_fails" is accurate. It was created in excel.
        dataIn_sub <- dataIn %>%
                select(progress_state,
                       global_pass,
                       max_training_sess_pass,
                       min_perf_pass,
                       global_pass_incl_phase_2_fails)
        
        dataIn_sub$compare_column <- dataIn_sub$global_pass == 1 | ( 
                dataIn_sub$progress_state == 'qc_failed_phase_2' & 
                        (dataIn_sub$max_training_sess_pass == 0 | 
                                 dataIn_sub$min_perf_pass == 0)
        )
        
        # Is there any mismatch? If so, make an error
        if (!identical(as.numeric(dataIn_sub$compare_column),
                       dataIn_sub$global_pass_incl_phase_2_fails)){
                stop('The column global_pass_incl_phase_2_fails is inaccurate')
        }
        
        # 2.
        
}