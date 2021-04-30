%% Draft script to analyze study 2 data

% This script will load the already analyzed files from the
% analyze_dependent_variables.m:
% results_table_all_ptp_analyzed.mat and
% long_form_data_all_ptp_analyzed.mat, and will add some additional columns
% such as global_pass_incl_phase_2_fail, columns from prolific meta data,
% etc.

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

%% For the long form data, add the current concept, curr arr, prompt_point_idx
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

% Now, prompt point idx

% - reshape...
item_point_idxs_mat = reshape(cell2mat(long_form_data_all_ptp.item_point_idxs),2,height(long_form_data_all_ptp))';

idx_prompt_item_idx_1 = long_form_data_all_ptp.prompt_item_idx == 1;
idx_prompt_item_idx_2 = long_form_data_all_ptp.prompt_item_idx == 2;

long_form_data_all_ptp.prompt_point_idx(idx_prompt_item_idx_1) = item_point_idxs_mat(idx_prompt_item_idx_1,1);
long_form_data_all_ptp.prompt_point_idx(idx_prompt_item_idx_2) = item_point_idxs_mat(idx_prompt_item_idx_2,2);

%% Add the prolific metadata columns

% Load the prolific metadata
prolific_metadata = readtable(fullfile(home,'data','prolific_meta_data','united_meta_data.xlsx'));

% Load the file with the mappings between entered and real prolific IDs....
entered_real_prolific_ids = load(fullfile(home,'results','analysis',...
    'other','entered_vs_real_prolific_ids.mat'));
entered_real_prolific_ids = entered_real_prolific_ids.entered_real_prolific_ids;

prolific_metadata_rel = table;
prolific_metadata_rel = prolific_metadata(:,{'participant_id',...
    'status','started_datetime','completed_date_time','time_taken',...
    'time_taken_min','age','entered_code','Sex'});

for iptp = 1:height(results_table_all_ptp)
    
    entered_id = results_table_all_ptp.ptp{iptp};
    
    if strcmp(entered_id,'test')
        
        results_table_all_ptp.participant_id{iptp} = '';
        
        continue
    end
    
    entered_id_prolific_idx = ...
        find(strcmp(entered_real_prolific_ids.entered_prolific_id,entered_id));
    
    if isempty(entered_id_prolific_idx)
        
        assert(strcmp(results_table_all_ptp.progress_state{iptp},'condition not assigned'));
        
        results_table_all_ptp.participant_id{iptp} = '';
        
        continue
        
    end
    
    assert(length(entered_id_prolific_idx) == 1);
    
    real_id = entered_real_prolific_ids.real_prolific_id{entered_id_prolific_idx};
    
    real_id_prolific_idx = find(strcmp(prolific_metadata_rel.participant_id,real_id));
    
    % Add the real id
    results_table_all_ptp.participant_id{iptp} = real_id;
    
    
end
%     prolific_rows = prolific_metadata_rel(real_id_prolific_idx,:);
%     
%     % Go through each column, and add to the results table
%     for iCol = 1:width(prolific_rows)
%         
%         curr_col = prolific_rows.Properties.VariableNames{iCol};
%         
%         results_table_all_ptp.([curr_col]){iptp} = prolific_rows{1,iCol};
%         
%         
%     end
%     
% end
    
    
