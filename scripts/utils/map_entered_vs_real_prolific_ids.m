%% Description

% Date: 20 April, 2021

% This script will take the real prolific IDs from the metadata file
% downloaded from prolific.co, and for each it will try to find a
% corresponding participant in the list of entered prolific IDs.

% The reason for making this script was that some participants messed up
% how they copy-pasted their prolific IDs, so their IDs in my raw data dont
% match their real prolific IDs. This created some issues in the pipeline,
% when trying to integrated prolific metadata with preprocessed and
% analysed data.

% This script is written after all the data was preprocessed, so its a
% post-hoc script. A new script should be written if we are going to get
% more data, such that preprocessing pipeline includes this mapping step.

%% General setup
clear; close all;
dbstop if error;

% Check you're in the right directory 
home = pwd;
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end
addpath(genpath(home));

%% Load the prolific meta data and the entered Prolific IDs from our preprocessed file

prolific_meta_data = readtable(fullfile(home,'data','prolific_meta_data',...
    'united_meta_data.xlsx'));

load(fullfile(home,'results','analysis','results_table_all_ptp.mat'));

%% Start the loop for each prolific meta data participant

% Insert a new columns 
other_file_row_idx = nan(height(prolific_meta_data),1);
prolific_meta_data = addvars(prolific_meta_data,other_file_row_idx,'Before',...
    'started_datetime');

other_file_row_idx = nan(height(results_table_all_ptp),1);
results_table_all_ptp = addvars(results_table_all_ptp,other_file_row_idx,'Before',...
    'maxPhase');

for ipptp = 1:height(prolific_meta_data)
    
    curr_ptp = prolific_meta_data.participant_id{ipptp};
    
    logical_strfind = strfind(results_table_all_ptp.ptp,curr_ptp);
    logical_strfind = cellfun(@isempty,logical_strfind,'UniformOutput',false);
    logical_strfind = cell2mat(logical_strfind);
    ind_match = find(logical_strfind == 0);

    if length(ind_match) == 2
       
        % So the participant probably attempted the experiment at first,
        % didn't work, and then did it again.
        
        matching_rows = results_table_all_ptp(ind_match,:);
        
        idx_row_we_want = find(~strcmp(matching_rows.progress_state,...
            'condition not assigned'));
        
        ind_match = ind_match(idx_row_we_want);
        
    end
    
    if ~isempty(ind_match)
        results_table_all_ptp.other_file_row_idx(ind_match) = ipptp;
        prolific_meta_data.other_file_row_idx(ipptp) = ind_match;
    end
    
end

%% Deal with any miscopied names

% So at this point, we can check the results_table_all_ptp and look for
% those that don't have a corresponding id in prolific meta data AND their
% status is NOT 'condition not assigned'. I find only 3 such participants,
% and for them we have to look for the rows in prolific meta data.
% Whatever remains in prolific metadata without matches, must have not
% submitted any data to any components of the experiment.

% Find the ones without a match in results table
idx_nan = find(isnan(results_table_all_ptp.other_file_row_idx));

% Find the ones without the progress state of 'condition not assigned'
idx_NOT_cond_not_ass = find(~strcmp(results_table_all_ptp.progress_state,...
    'condition not assigned'));

% also exclude the 'test' row, which was just me debugging the script
idx_not_test = find(~strcmp(results_table_all_ptp.ptp,'test'));

idx_nan_NOT_cond_not_ass_not_test = intersect(...
    intersect(idx_nan,idx_NOT_cond_not_ass),idx_not_test);

% For these 3 participants, find corresponding prolific meta data rows
for iMiss = 1:length(idx_nan_NOT_cond_not_ass_not_test)
    
    iRow = idx_nan_NOT_cond_not_ass_not_test(iMiss);
    
    fprintf(['Searching for ' results_table_all_ptp.ptp{iRow} ...
        ' in prolific meta data...\n']);
    
    % Take first 8 characters
    str_to_search = results_table_all_ptp.ptp{iRow}(1:8);
    
    % Find this in prolific
    logical_strfind = strfind(prolific_meta_data.participant_id,str_to_search);
    logical_strfind = cellfun(@isempty,logical_strfind,'UniformOutput',false);
    logical_strfind = cell2mat(logical_strfind);
    ind_match = find(logical_strfind == 0);    
    
    if ~isempty(ind_match)
        
        fprintf(['Found! Row ' int2str(ind_match) ...
            '. Full ID: ' prolific_meta_data.participant_id{ind_match} '\n']);
    
        results_table_all_ptp.other_file_row_idx(iRow) = ind_match;
        
        prolific_meta_data.other_file_row_idx(ind_match) = iRow;
        
    else
        fprintf('Failed...\n');
    end
end

%% Save everything in a separate file
idx_rows = find(~isnan(results_table_all_ptp.other_file_row_idx));

entered_real_prolific_ids = results_table_all_ptp(idx_rows,...
    {'ptp','other_file_row_idx'});

real_prolific_id = prolific_meta_data.participant_id(...
    entered_real_prolific_ids.other_file_row_idx);

entered_real_prolific_ids = addvars(entered_real_prolific_ids,...
    real_prolific_id);

entered_real_prolific_ids.other_file_row_idx = [];

entered_real_prolific_ids = renamevars(entered_real_prolific_ids,...
    'ptp','entered_prolific_id');

%% Save this file
if ~exist(fullfile(home,'results','analysis','other'))
    mkdir(fullfile(home,'results','analysis','other'));
end

writetable(entered_real_prolific_ids,fullfile(home,'results','analysis',...
    'other','entered_vs_real_prolific_ids.csv'));

save(fullfile(home,'results','analysis','other',...
    'entered_vs_real_prolific_ids.mat'),'entered_real_prolific_ids');