# Description

# This function takes in one of the results tables, and deletes data from those 
# participants who did not pass the QC criteria

# Potentially, add a condition that the participant did not "RETURN" 
# their prolific submission

get_only_qc_pass = function(dataIn){
        
        dataIn <- subset(dataIn,
                                 data_submitted                 == 1 & 
                                 global_pass_incl_phase_2_fails == 1 & 
                                 debrief_qc_pass                == 1 & 
                                 fb_int_qc_pass                 == 1 & 
                                 phase_1_rt_qc_pass             == 1 & 
                                 phase_2_rt_qc_pass             == 1 &
                                 basic_data_checks_pass         == 1)

        return(dataIn)
}