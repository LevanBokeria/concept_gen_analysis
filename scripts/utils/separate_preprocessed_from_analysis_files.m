%% Description

% Will load the results tables, delete those columns that
% were added as analysis columns, such as learning rate exponent etc, so
% all columns that were not part of the preprocessing script pipeline.
% The idea is to just keep these preprocessed files as untouched, and have
% a separate script that then loads those, adds any new columns, and saves
% as a different file. 

%% General setup
clear; clc;


% Check you're in the right directory 
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

% Flags

saveFiles = 1;

%% Load the data for all ptps:
[~,results_table_all_ptp,...
    ~, ~] = ...
    load_all_participant_files(home);

%% Delete the extra columns

cols_to_delete = {'phase_1_learning_rate_exp','phase_2_learning_rate_exp',...
    'phase_1_learning_rate_int','phase_2_learning_rate_int',...
    'phase_1_learning_rate_sse','phase_2_learning_rate_sse',...
    'phase_2_min_phase_1_learning_rate_exp',...
    'phase_2_min_phase_1_ses_1_2_perf'};

results_table_all_ptp = removevars(results_table_all_ptp,cols_to_delete);

%% Save the new variable
if saveFiles
    save(fullfile(home,'results','analysis','results_table_all_ptp.mat'),...
        'results_table_all_ptp')
end














