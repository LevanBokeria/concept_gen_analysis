function [long_form_data,results_table] = ...
    basic_data_checks(long_form_data,results_table)

%% Description
% Check ptp data for basic issues, such as n trials, n sessions, etc.
% Takes in the long_form data.

%% Basic setup
dbstop if error;

if (nargin < 1)
    
    clear; clc;
    
    % Get the current directory
    home = pwd;
    
    [~,name,~] = fileparts(home);
    
    if ~strcmp('concept_gen_analysis',name)
        error('please change working directory to ./con_learn/concept_gen_analysis/');
    end
    
    % Load the files
    [long_form_data,results_table,...
        ~, ~] = ...
        load_all_participant_files(home);
    
end

%% Start the data checks

% Add a column saying whether the basic data checks were violated or not
long_form_data.basic_data_checks_pass = ...
    true(height(long_form_data),1);
results_table.basic_data_checks_pass = ...
    true(height(results_table),1);

% - n trials for each session

G = groupsummary(long_form_data,{'ptp','phase','session'});

% Did anyone have trials more than 42?
idx_too_many_trials = find(G.GroupCount > 42);

if ~isempty(idx_too_many_trials)
    
    fprintf(['Found ' int2str(numel(idx_too_many_trials))...
        ' participants with ' ...
        'too many trials...\n']);
    
    % Get their participant names
    ptp_too_many_trials = G.ptp(idx_too_many_trials);
    
    for iP = 1:numel(ptp_too_many_trials)
        
        curr_ptp = ptp_too_many_trials{iP};
        
        % Find this in the long form
        idx_match_long = find(strcmp(long_form_data.ptp,curr_ptp));
        
        long_form_data.basic_data_checks_pass(idx_match_long) = false;
        
        % Find in the results table
        idx_match_results = find(strcmp(results_table.ptp,curr_ptp));
        results_table.basic_data_checks_pass(idx_match_results) = false;
        
    end
    
end % if ~isempty idx


end