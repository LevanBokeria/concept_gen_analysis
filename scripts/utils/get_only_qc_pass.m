function data_out = get_only_qc_pass(data_in)

%% Description

% This function will just filter the given dataset to give you only data
% with qc passed participants.

% It will also include those participants that failed QC in phase 2 but
% because of not reaching performance thresholds.

%% If no arguments
if (nargin < 1)
    
    % Get the current directory
    home = pwd;
    
    [~,name,~] = fileparts(home);
    
    if ~strcmp('concept_gen_analysis',name)
        error('please change working directory to ./con_learn/concept_gen_analysis/');
    end
    
    % Load the files
    [long_form_data_all_ptp,results_table_all_ptp,...
        debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
        load_all_participant_files(home);
    
    data_in = results_table_all_ptp;
    
end

%% Clean the data

% Get the ones that had global pass == 1, OR if failed at phase 2 because
% of performance, include them too.
idx_global_pass = data_in.global_pass == 1;
idx_phase_2_perf_fail = strcmp(data_in.progress_state,'qc_failed_phase_2') & ...
    (data_in.min_perf_pass == 0 | data_in.max_training_sess_pass == 0);
idx_global_pass_or_phase_2_perf_fail = idx_global_pass | idx_phase_2_perf_fail;

idx_qc_pass_up_to_now = idx_global_pass_or_phase_2_perf_fail & ...
    data_in.debrief_qc_pass == 1 & ...
    data_in.fb_int_qc_pass == 1 & ...
    data_in.data_submitted == 1;

data_out = data_in(idx_qc_pass_up_to_now,:);

all_ptps = unique(data_out.ptp);

end