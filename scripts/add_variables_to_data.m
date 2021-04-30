%% Draft script to analyze study 2 data

% This script will load the already analyzed files
% results_table_all_ptp_analyzed.mat and
% long_form_data_all_ptp_analyzed.mat, and will add some additional columns
% such as global_pass_incl_phase_2_fail, columns from prolific meta data,
% time taken in minutes, etc etc.

%% Global parameters
clear; close all;

dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

addpath(genpath(home));

%% Data parameters 
saveFiles = 0;

%% Load the data for all ptps:
load(fullfile(home,'results','analysis','results_table_all_ptp_analyzed.mat'));
load(fullfile(home,'results','analysis','long_form_data_all_ptp_analyzed.mat'));

nPtp = height(results_table_all_ptp);

%% Add the "global pass including phase 2 fail" column, to both files
% This column will check if a participant has global_pass=0 but they failed
% in phase 2 because of performance checks. Such participants should still
% be analyzed even thought they didn't reach the criterion in phase 2. 

% Add these columns to the data
long_form_data_all_ptp.global_pass_incl_phase_2_fails = long_form_data_all_ptp.global_pass;
results_table_all_ptp.global_pass_incl_phase_2_fails  = results_table_all_ptp.global_pass;

% Find indices
idx_pass_long_form = long_form_data_all_ptp.global_pass == 0 & ...
    strcmp(long_form_data_all_ptp.progress_state,'qc_failed_phase_2') & ...
    (long_form_data_all_ptp.min_perf_pass == 0 | ...
    long_form_data_all_ptp.max_training_sess_pass == 0);

idx_pass_results_table = results_table_all_ptp.global_pass == 0 & ...
    strcmp(results_table_all_ptp.progress_state,'qc_failed_phase_2') & ...
    (results_table_all_ptp.min_perf_pass == 0 | ...
    results_table_all_ptp.max_training_sess_pass == 0);

% Change the values
long_form_data_all_ptp.global_pass_incl_phase_2_fails(idx_pass_long_form) = 1;
results_table_all_ptp.global_pass_incl_phase_2_fails(idx_pass_results_table) = 1;

%% For the long form data, add the current concept, curr arr, prompt_point_idx and time elapsed in mins
long_form_data_all_ptp.current_concept     = cell(height(long_form_data_all_ptp),1);
long_form_data_all_ptp.arr_phase_1_name    = zeros(height(long_form_data_all_ptp),1);
long_form_data_all_ptp.arr_phase_2_name    = zeros(height(long_form_data_all_ptp),1);
long_form_data_all_ptp.current_arrangement = zeros(height(long_form_data_all_ptp),1);
long_form_data_all_ptp.prompt_point_idx    = zeros(height(long_form_data_all_ptp),1);
long_form_data_all_ptp.time_elapsed_min    = zeros(height(long_form_data_all_ptp),1);


% Change arr_phase_1_name and phase 2

% Phase 1 arr 1 idxs
arr_phase_1_mat = reshape(cell2mat(long_form_data_all_ptp.arr_phase_1),3,height(long_form_data_all_ptp))';
arr_phase_2_mat = reshape(cell2mat(long_form_data_all_ptp.arr_phase_2),3,height(long_form_data_all_ptp))';

% Find all the indices to substitute
arr_1_phase_1_idx = arr_phase_1_mat(:,1) == 1;
arr_2_phase_1_idx = arr_phase_1_mat(:,1) == 14;
arr_3_phase_1_idx = arr_phase_1_mat(:,1) == 3;
arr_4_phase_1_idx = arr_phase_1_mat(:,1) == 15;

arr_1_phase_2_idx = arr_phase_2_mat(:,1) == 1;
arr_2_phase_2_idx = arr_phase_2_mat(:,1) == 14;
arr_3_phase_2_idx = arr_phase_2_mat(:,1) == 3;
arr_4_phase_2_idx = arr_phase_2_mat(:,1) == 15;

% Populate arrangement names 
long_form_data_all_ptp.arr_phase_1_name(arr_1_phase_1_idx) = 1;
long_form_data_all_ptp.arr_phase_1_name(arr_2_phase_1_idx) = 2;
long_form_data_all_ptp.arr_phase_1_name(arr_3_phase_1_idx) = 3;
long_form_data_all_ptp.arr_phase_1_name(arr_4_phase_1_idx) = 4;

long_form_data_all_ptp.arr_phase_2_name(arr_1_phase_2_idx) = 1;
long_form_data_all_ptp.arr_phase_2_name(arr_2_phase_2_idx) = 2;
long_form_data_all_ptp.arr_phase_2_name(arr_3_phase_2_idx) = 3;
long_form_data_all_ptp.arr_phase_2_name(arr_4_phase_2_idx) = 4;

% Current arrangement
% - find all phase 1 idxs
idx_phase_1 = long_form_data_all_ptp.phase == 1;
idx_phase_2 = long_form_data_all_ptp.phase == 2;

long_form_data_all_ptp.current_arrangement(idx_phase_1) = long_form_data_all_ptp.arr_phase_1_name(idx_phase_1);
long_form_data_all_ptp.current_arrangement(idx_phase_2) = long_form_data_all_ptp.arr_phase_1_name(idx_phase_2);


% Do the same for current concept
long_form_data_all_ptp.current_concept(idx_phase_1) = long_form_data_all_ptp.concept_phase_1(idx_phase_1);
long_form_data_all_ptp.current_concept(idx_phase_2) = long_form_data_all_ptp.concept_phase_2(idx_phase_2);

