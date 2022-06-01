# Preprocessing pipeline in Matlab:

1. In the preprocessing folder, start by using the appropriate "preprocessing_wave....m" files, which will create preprocessed folders and files for all the participants specified.

2. Go through the wave debriefing, score whether each participant pass/failed the QC based on debriefing feedback.

3. Run the debrief_and_int_feedback_analysis_updating.m script, to take the variables from excel about pass/fail in debriefs, and add them to the wave_debriefing files themselves.

4. Run the wrapper_add_recalc_how_many.m script. This will add the wave data to the combined data of all the participants, then it will recalculate outliers based on RT, then it will tell you how many more participants you need to complete data acquisition for a given batch of the wave. You use this to open appropriate number of slots in prolific.

# Analysis pipeline in Matlab:

1. analyze_dependent_variables.m this will fit the inverse exponential to get the learning rate for each participant. It will also calculate phase 2 minus phase 1 scores for each participant for both of the dependent variables - learning rate and ses1-2 performance, and add those to the results_table_all_ptp table. 

2. add_variables.to_data.m this will load results_table_all_ptp_analyzed.mat and long_form_data_all_ptp_analyzed.mat. It will add the following columns to these tables: 
- global_pass_including_phase_2_fail (both tables).
- current_concept, current_arrangement, arr_phase_1_name, arr_phase_2_name, prompt_point_idx, time_elapsed_min (long_form table).
- columns from the prolific_meta_data file (results_table).

# Analysis continued in R: