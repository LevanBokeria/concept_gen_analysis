function wrapper_add_recalc_how_many(batch_number,wave_name,...
    overwriteData,saveCombinedData)

% Description:

% This is a wrapper script, to quickly call the following functions after
% having analyzed specific wave data:

% - add_to_all_ptp.m
% - recalc_rt_qc_pass_func.m
% - how_many_new_ptp_needed.m

%% Global variables

if (nargin < 1)
    batch_number    = 'batch_3';
    wave_name       = 'wave_17_ds';
    overwriteData    = 0;
    saveCombinedData = 1;
end

% Get the current directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

batch_wave_name = [batch_number '_' wave_name];

%% Add to all ptps
add_to_all_ptp(batch_number,wave_name,overwriteData,saveCombinedData);

%% Recalculate RT QC pass checks
[long_form_data_all_ptp,results_table_all_ptp] = ...
    recalc_rt_qc_pass_func(saveCombinedData);

%% How many participants are needed?
[long_form_data_all_ptp,results_table_all_ptp,...
    G,IID,...
    results_table_qc_pass,table_for_groups] = ...
    how_many_new_ptp_needed;

%% Save all the excel files

% % Load the files
% [long_form_data_all_ptp,results_table_all_ptp,...
%     debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
%     load_all_participant_files(home);
% % Load the qc pass analyzed stuff



end