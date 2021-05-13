%% Draft script to analyze study 2 data

% This script will load the already analyzed files from the
% analyze_dependent_variables.m:
% results_table_all_ptp_analyzed.mat and
% long_form_data_all_ptp_analyzed.mat, and will add some additional columns
% such as global_pass_incl_phase_2_fail, columns from prolific meta data,
% etc.

% It will also call a script that checks for some basic data checks, such
% as that number of trials, sessions, phases, etc was correct for each
% participant.

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
saveLongFormMat = 0;
saveLongFormCSV = 0;
saveResultsTableMat = 1;
saveResultsTableCSV = 1;

%% Load the data for all ptps:
load(fullfile(home,'results','analysis','results_table_all_ptp_analyzed.mat'));
load(fullfile(home,'results','analysis','long_form_data_all_ptp_analyzed.mat'));

nPtp = height(results_table_all_ptp);

%% Add the "global pass including phase 2 fail" column, to both files
% This column will check if a participant has global_pass=0 but they failed
% in phase 2 because of performance checks. Such participants should still
% be analyzed even thought they didn't reach the criterion in phase 2. 

fprintf('Adding the global_pass_incl_phase_2_fails columns...\n');

% % If the column already exists, alert me and skip this
% if ~isempty(find(strcmp(results_table_all_ptp.Properties.VariableNames,...
%         'global_pass_incl_phase_2_fails')))
%     fprintf('results table already has the global_pass_incl_phase_2_fails column!!! Skipping this step...\n');
% else
    
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
    long_form_data_all_ptp.global_pass_incl_phase_2_fails(idx_pass_long_form)    = 1;
    results_table_all_ptp.global_pass_incl_phase_2_fails(idx_pass_results_table) = 1;

% end
%% For the long form data, add the current concept, curr arr, prompt_point_idx

fprintf(['Adding current_concept, current_arrangement, ',...
    'prompr_point_idx, time_elapsed_min, arr_phase_1_name and ', ...
    'arr_phase_2_name columns to the long form data table...\n']);

% Check if these already exist
% if ~isempty(find(strcmp(long_form_data_all_ptp.Properties.VariableNames,...
%         'experiment')))
%     
%     fprintf('long_form_data file already has a current_concept column. Probably others too. Skipping this step...\n');
% else
    long_form_data_all_ptp.current_concept     = cell(height(long_form_data_all_ptp),1);
    long_form_data_all_ptp.arr_phase_1_name    = zeros(height(long_form_data_all_ptp),1);
    long_form_data_all_ptp.arr_phase_2_name    = zeros(height(long_form_data_all_ptp),1);
    long_form_data_all_ptp.current_arrangement = zeros(height(long_form_data_all_ptp),1);
    long_form_data_all_ptp.prompt_point_idx    = zeros(height(long_form_data_all_ptp),1);
    long_form_data_all_ptp.experiment          = zeros(height(long_form_data_all_ptp),1);
    
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
    long_form_data_all_ptp.current_arrangement(idx_phase_2) = long_form_data_all_ptp.arr_phase_2_name(idx_phase_2);
    
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

    % Now the experiment
    long_form_data_all_ptp.experiment(arr_1_phase_1_idx | arr_2_phase_1_idx) = ...
        1;
    long_form_data_all_ptp.experiment(arr_3_phase_1_idx | arr_4_phase_1_idx) = ...
        2;    
% end

%% Add the arrangement names and experiment to the results table
results_table_all_ptp.arr_phase_1_name    = zeros(height(results_table_all_ptp),1);
results_table_all_ptp.arr_phase_2_name    = zeros(height(results_table_all_ptp),1);
results_table_all_ptp.experiment          = zeros(height(results_table_all_ptp),1);

for i = 1:height(results_table_all_ptp)
    
    % Give names to arrangements, and to the experiment
    if results_table_all_ptp.arr_phase_1{i} == [1;8;15]
        results_table_all_ptp.arr_phase_1_name(i) = 1;
        
        if results_table_all_ptp.congruency(i) == 1
            results_table_all_ptp.arr_phase_2_name(i) = 1;
        else
            results_table_all_ptp.arr_phase_2_name(i) = 2;
        end
        
        results_table_all_ptp.experiment(i)       = 1;
        
    elseif results_table_all_ptp.arr_phase_1{i} == [14;9;4]
        results_table_all_ptp.arr_phase_1_name(i) = 2;

        if results_table_all_ptp.congruency(i) == 1
            results_table_all_ptp.arr_phase_2_name(i) = 2;
        else
            results_table_all_ptp.arr_phase_2_name(i) = 1;
        end        
        results_table_all_ptp.experiment(i)       = 1;        
        
    elseif results_table_all_ptp.arr_phase_1{i} == [3;8;14]
        results_table_all_ptp.arr_phase_1_name(i) = 3;

        if results_table_all_ptp.congruency(i) == 1
            results_table_all_ptp.arr_phase_2_name(i) = 3;
        else
            results_table_all_ptp.arr_phase_2_name(i) = 4;
        end        
        results_table_all_ptp.experiment(i)       = 2;        
        
    elseif results_table_all_ptp.arr_phase_1{i} == [15;5;12]
        results_table_all_ptp.arr_phase_1_name(i) = 4;

        if results_table_all_ptp.congruency(i) == 1
            results_table_all_ptp.arr_phase_2_name(i) = 4;
        else
            results_table_all_ptp.arr_phase_2_name(i) = 3;
        end        
        results_table_all_ptp.experiment(i)       = 2;        
        
    elseif isnan(results_table_all_ptp.arr_phase_1{i})
        results_table_all_ptp.arr_phase_1_name(i) = NaN;
        results_table_all_ptp.arr_phase_2_name(i) = NaN;
        
        results_table_all_ptp.experiment(i)       = NaN;
        
    else
        i
    end
    
end
%% Add the prolific metadata columns to the results table

fprintf('Adding prolific meta data to the results table... \n');


% Check if prolific data is already combined
if ~isempty(find(strcmp(results_table_all_ptp.Properties.VariableNames,...
        'status')))
    fprintf('results table seems to be already combined with prolific meta data. Skipping... \n');
else
    
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
    
    % We're going to use the join function to combine the data in the two
    % tables. For that, we need a shared "key" column. So in the results_table
    % create the participant_id column, the entered prolific IDs.
    
    % For those rows that don't have a corresponding entered ID, remove them.
    
    idx_to_remove = [];
    
    for iptp = 1:height(results_table_all_ptp)
        
        entered_id = results_table_all_ptp.ptp{iptp};
        
        if strcmp(entered_id,'test')
            
            results_table_all_ptp.participant_id{iptp} = entered_id;
            
            idx_to_remove = [idx_to_remove iptp];
            
            % Add this as an empty row to the prolific meta data
            prolific_metadata_rel.participant_id{height(prolific_metadata_rel)...
                +1} = entered_id;
            
            continue
        end
        
        entered_id_prolific_idx = ...
            find(strcmp(entered_real_prolific_ids.entered_prolific_id,entered_id));
        
        if isempty(entered_id_prolific_idx)
            
            assert(strcmp(results_table_all_ptp.progress_state{iptp},'condition not assigned'));
            
            fprintf(['Participant ' entered_id ' not found in the prolific meta data.\n']);
            
            results_table_all_ptp.participant_id{iptp} = entered_id;
            
            idx_to_remove = [idx_to_remove iptp];
            
            % Add this as an empty row to the prolific meta data
            prolific_metadata_rel.participant_id{height(prolific_metadata_rel)...
                +1} = entered_id;
            
            continue
            
        end
        
        assert(length(entered_id_prolific_idx) == 1);
        
        real_id = entered_real_prolific_ids.real_prolific_id{entered_id_prolific_idx};
        
        real_id_prolific_idx = find(strcmp(prolific_metadata_rel.participant_id,real_id));
        
        % Add the real id
        results_table_all_ptp.participant_id{iptp} = real_id;
        
    end
    
    % Join the tables
    results_table_all_ptp = join(results_table_all_ptp,prolific_metadata_rel,...
        'Keys','participant_id');
    
    % remove the participant ID column
    results_table_all_ptp.participant_id = [];

end

%% Do the basic data checks 
fprintf('Doing the basic data checks now...\n');
[long_form_data_all_ptp,results_table_all_ptp] = ...
    basic_data_checks(long_form_data_all_ptp,results_table_all_ptp);

%% Save everything

if saveLongFormMat
    
    fprintf('Saving long form mat file...\n');

    save(fullfile(home,'results','analysis',...
        'long_form_data_all_ptp_analyzed.mat'),'long_form_data_all_ptp');
end

if saveResultsTableMat
    
    fprintf('Saving results table mat file...\n');
    
    save(fullfile(home,'results','analysis',...
        'results_table_all_ptp_analyzed.mat'),'results_table_all_ptp');
end

if saveLongFormCSV
    
    fprintf('Saving long form csv file...\n');

    save_table_for_excel(long_form_data_all_ptp,fullfile(home,'results',...
        'analysis','long_form_data_all_ptp_analyzed.csv'),1);    
end

if saveResultsTableCSV
    
    fprintf('Saving results table csv file...\n');

    save_table_for_excel(results_table_all_ptp,fullfile(home,'results',...
        'analysis','results_table_all_ptp_analyzed.csv'),1);

end