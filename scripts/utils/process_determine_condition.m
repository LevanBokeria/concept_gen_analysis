%% Description

% This script will analyze just the determine_condition module output, and
% tell you what condition that participant had.
% This was used to quickly tell the conditions of all the people that were
% ran up to a certain point.

% Then, the script will read an excel file that is the output of prolific,
% telling who failed QC, who returned the submission, and who is ongoing,
% or finished successfully.

% Then, the script will take only those participants that are ongoing or
% passed, and tabulate how many have what consitions, so we can determine
% who else to run. 

%% Global variables
clear; clc;

% Get the current directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

%% Read the filename and determine the number of participants
jatos_file_name    = 'jatos_results_20210222203331';
prolific_file_name = 'prolific_export_600eca365ae8f07c9874d9cd';

%% First, read the prolific data table
prolific_data = readtable(fullfile(home,'data','sandbox',...
    'determine_condition_files',[prolific_file_name '.csv']));

%% Now, read the jatos determine condition file, and create a results table
rawtext = fileread(fullfile(home,'data','sandbox','determine_condition_files',...
    [jatos_file_name '.txt']));

% Split at the first component
[ptp_split,~] = split(rawtext,'[conditions_start---');

nPtp = numel(ptp_split) - 1;

results_table = table;

for iPtp = 1:nPtp
    
    % Take this data
    curr_cell = ptp_split{iPtp+1};
    
    % Get the data submission section
    splitter_end = '---conditions_end]';
    [curr_cell,matches] = split(curr_cell,splitter_end);
    
    decoded = jsondecode(curr_cell{1});
    
    %% Get the variables to record
    
    iCongruency             = decoded.inputData.congruency;
    iConcept_phase_1        = decoded.inputData.concepts.phase_1.concept_space;
    iConcept_phase_2        = decoded.inputData.concepts.phase_2.concept_space;
    iArr_phase_1            = decoded.inputData.basic_parameters.targetPoints.phase_1;
    iArr_phase_2            = decoded.inputData.basic_parameters.targetPoints.phase_2;
    
    targetNamesUsed_phase_1 = decoded.inputData.basic_parameters.targetNamesUsed.phase_1;
    targetNamesUsed_phase_2 = decoded.inputData.basic_parameters.targetNamesUsed.phase_2;
    
    iProgressState          = decoded.progress_state;
    
    %% Put all of this in a table
    curr_ptp                                    = decoded.prolific_ID;
    
    % Some ptps entered their IDs incorrectly
    if strcmp(curr_ptp,'5ffc4646cecb15016da849e5ffc4646cecb15016da849e5')
        curr_ptp = '5ffc4646cecb15016da849e5';
    elseif strcmp(curr_ptp,'59b70c4e7547b100012d6d4259b70c4e7547b100012d6d42')
        curr_ptp = '59b70c4e7547b100012d6d42';
    elseif strcmp(curr_ptp,'5fc3786339c5f2192e8c5fc3786339c5f2192e8c129c')
        curr_ptp = '5fc3786339c5f2192e8c129c';
    end
    
    results_table.ptp{iPtp}                     = curr_ptp;
    results_table.progress_state{iPtp}          = iProgressState;
      
    % Get the details of the condition
    results_table.congruency(iPtp)              = iCongruency;
    results_table.concept_phase_1{iPtp}         = iConcept_phase_1;
    results_table.concept_phase_2{iPtp}         = iConcept_phase_2;
    results_table.arr_phase_1{iPtp}             = iArr_phase_1;
    results_table.arr_phase_2{iPtp}             = iArr_phase_2;

    %% Cross check this participant in prolific data file
    
    % Find the row
    find_row = cellfun(@strfind,prolific_data.participant_id,...
        repmat({curr_ptp},height(prolific_data),1),'UniformOutput',false);
    
    find_empty = cellfun(@isempty,find_row);
    
    idx_row = find(find_empty == 0);
    
    if isempty(idx_row)
        curr_status = 'not_in_prolific_table';       
        entered_code = NaN;
    else
        assert(numel(idx_row) == 1);
        curr_status = prolific_data.status{idx_row};
        
        entered_code = prolific_data.entered_code{idx_row};
    end
    
    results_table.prolific_status{iPtp} = curr_status;
    results_table.entered_code{iPtp}    = entered_code;

end

%% Filter out anything that is returned or code is QC fail

find_ret = cellfun(@strfind,results_table.prolific_status,...
    repmat({'RETURNED'},nPtp,1),'UniformOutput',false);
find_ret_idx = find(cellfun(@isempty,find_ret));

% Find qc fails
find_qc_fail = cellfun(@strfind,results_table.entered_code,...
    repmat({'QC_fail'},nPtp,1),'UniformOutput',false);
find_qc_fail_idx = find(cellfun(@isempty,find_qc_fail));

find_ret_qc_fail_idx = intersect(find_ret_idx,find_qc_fail_idx);

% Find timeouts
find_timeout = cellfun(@strfind,results_table.prolific_status,...
    repmat({'TIMED-OUT'},nPtp,1),'UniformOutput',false);
find_timeout_idx = find(cellfun(@isempty,find_timeout));

find_ret_qc_fail_timeout_idx = intersect(find_ret_qc_fail_idx,find_timeout_idx);

% Find practice fail
find_practice_fail = cellfun(@strfind,results_table.entered_code,...
    repmat({'PRACTICE_FAILED'},nPtp,1),'UniformOutput',false);
find_practice_fail_idx = find(cellfun(@isempty,find_practice_fail));

find_ret_qc_fail_timeout_practice_fail_idx = ...
    intersect(find_ret_qc_fail_timeout_idx,find_practice_fail_idx);

results_table_clean = results_table(find_ret_qc_fail_timeout_practice_fail_idx,:);

%% Now, tabulate the conditions


% table to analyze
table_for_groups = results_table_clean(:,{'congruency','concept_phase_1',...
    'concept_phase_2','arr_phase_1','arr_phase_2'});

% Turn the arrangements to a string, so findgroups can classify them
table_for_groups.arr_phase_1 = cellfun(@transpose,...
    table_for_groups.arr_phase_1,'UniformOutput',0);
table_for_groups.arr_phase_1 = cellfun(@num2str,...
    table_for_groups.arr_phase_1,'UniformOutput',0);
table_for_groups.arr_phase_2 = cellfun(@transpose,...
    table_for_groups.arr_phase_2,'UniformOutput',0);
table_for_groups.arr_phase_2 = cellfun(@num2str,...
    table_for_groups.arr_phase_2,'UniformOutput',0);

[G,IID] = findgroups(table_for_groups);

% Add the counts to the table
for iRow = 1:height(IID)
    IID.count(iRow) = sum(G == iRow);
end











