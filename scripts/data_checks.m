%% Draft script to analyze study 2 data

% Check ptp data for basic issues

clear; close all;

dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

addpath(genpath(home));

task = 'pilots';

% Suppress warnings
warning('off','MATLAB:table:RowsAddedExistingVars')

participants = table;
% 
participants.names = {'pilot2','pilot3','pilot4','pilot5','pilot6',...
    'pilot7','pilot9','pilot10','pilot11','pilot15',...
    'pilot16','pilot17','pilot18','pilot19','pilot20',...
    'pilot21','pilot22','pilot23','pilot24','pilot26','pilot27'}';

% participants.names = {'pilot20'};

% participants.oldStyle     = [1,1,zeros(1,6)]';
participants.oldStyle     = [1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]';
% participants.oldStyle     = 0;

nPtp = height(participants);

savePerfFigure = 0;
averageAcrossPtp = 0;

% Flag to draw plots
drawPlots         = 0;
plotFMSEstimation = 0;
writeExcel        = 0;

% Create a big table to hold a long-format dataset of all ptp data
data_long = table;

for iPtp = 1:height(participants)
    
    % Reset the names variable
    names = {};
    
    curr_ptp = participants.names{iPtp}
    
    curr_style = participants.oldStyle(iPtp);
    
    %% Read and decode
    rawtext = fileread(fullfile(home,'data',task,['jatos_results_' curr_ptp '.txt']));

    % Split at the first component
    splitter_start = '[data_submission_start';
    [subject_split,matches] = split(rawtext,[splitter_start '---']);

    splitter_end = 'data_submission_end]';
    [split2,matches2] = split(subject_split{2},['---' splitter_end]);

    decoded = jsondecode(split2{1});

    participants.congruency(iPtp) = decoded.inputData.congruency;
    
    %% Get the time they took to read instructions
    participants.phase_1_instr_rt(iPtp) = decoded.outputData.instructions_results.phase_1.rt;
    try
        participants.phase_2_instr_rt(iPtp) = decoded.outputData.instructions_results.phase_2.rt;
    catch
        participants.phase_2_instr_rt(iPtp) = NaN;
    end
    try
        participants.phase_3_instr_rt(iPtp) = decoded.outputData.instructions_results.phase_3.rt;
    catch
        participants.phase_3_instr_rt(iPtp) = NaN;
    end
    
    
    if curr_style
        % Data is coded old way, separate struct for each phase

        % Show me quick relevant data
        qc_status = decoded.qc_status;

        % Running performance
        running_perf = decoded.inputData.running_perf;

        
        %% Start the phase loop
        phase_results = table;
        for iPhase = 1:3
            
            phase_tbl = decoded.outputData.(['phase_' int2str(iPhase) ...
                '_results']);
            
            nSessions = size(phase_tbl,1);
            
            for iSession = 1:nSessions
                
                ses_tbl = struct2table(phase_tbl(iSession,:));
                
                ses_tbl.ptp = cell(height(ses_tbl),1);
                [names{1:height(ses_tbl)}] = deal(curr_ptp);
                
                ses_tbl.ptp = names';
                
                % If missed, put NA there, for each column of data
                cols_to_touch = {'rt','key_press','correct'};
                
                for iCol = 1:numel(cols_to_touch)
                    curr_col = cols_to_touch{iCol};
                    
                    if isa(ses_tbl.(curr_col)(1),'cell') % so we have cells, theres a miss
                        idx_empty = find(cellfun(@isempty,ses_tbl.(curr_col)));
                        if ~isempty(idx_empty)
                            ses_tbl.(curr_col)(idx_empty) = {NaN};
                        end
                        ses_tbl.(curr_col) = cell2mat(ses_tbl.(curr_col));
                    end
                end
                
%                 % Change rt column to cell
%                 if isa(ses_tbl.rt(1),'double')
%                     ses_tbl.rt = ...
%                         num2cell(ses_tbl.rt);
%                 end
%                 
%                 if isa(ses_tbl.key_press(1),'double')
%                     ses_tbl.key_press = ...
%                         num2cell(ses_tbl.key_press);
%                 end               
% 
%                 if isa(ses_tbl.correct(1),'double')
%                     ses_tbl.correct = ...
%                         num2cell(ses_tbl.correct);
%                 end                               
                
                phase_results = [phase_results; ses_tbl];

            end
       
        end % iPhase
        
        % add the running performance and global trial columns, just to
        % concatenate with newer ptp data
        running_perf = cell(height(phase_results),1);
        phase_results = addvars(phase_results,running_perf,'Before','rt');
        
        global_trial = NaN(height(phase_results),1);
        phase_results = addvars(phase_results,global_trial,'After','internal_node_id');
        
    else
        % Data coded new way, in the row format 
        phase_results = struct2table(decoded.outputData.phase_results);
        
        % If missed, put NA there, for each column of data 
        cols_to_touch = {'rt','key_press','correct'};
        
        for iCol = 1:numel(cols_to_touch)
            curr_col = cols_to_touch{iCol};
            
            if isa(phase_results.(curr_col)(1),'cell') % so we have cells, theres a miss
                idx_empty = find(cellfun(@isempty,phase_results.(curr_col)));
                if ~isempty(idx_empty)
                    phase_results.(curr_col)(idx_empty) = {NaN};
                end
                phase_results.(curr_col) = cell2mat(phase_results.(curr_col));
            end
        end
        
        % If the practice column contains the word 'practice' change to -1
        % and turn the array to a double
        if isa(phase_results.session(1),'cell')
            idx_practice = strcmp(phase_results.session,'practice');
            if nnz(idx_practice) > 0
                phase_results.session(idx_practice) = {-1};
            end
            phase_results.session = cell2mat(phase_results.session);
        end
        
        % Add the participant name to the table
        phase_results.ptp = cell(height(phase_results),1);
        [names{1:height(phase_results)}] = deal(curr_ptp);
        
        phase_results.ptp = names';
        
    end % if curr_style

    %% Record the participant table
    participants.data{iPtp} = phase_results;
    
    data_long = [data_long; phase_results];
end % for iPtp

% Just a counter to write in the table
ctr = 1; 

%% Start the loop
for iPtp = 1:nPtp
    
    % Name of current participant
    curr_ptp = participants.names{iPtp}
    
    for iPhase = 1:2
        
        %% Get data for this ptp for this phase, minus practice trials
        idx_ptp      = strcmp(data_long.ptp,curr_ptp);
        idx_phase    = data_long.phase == iPhase;
        % Find practice trials
        idx_practice = data_long.session == -1;
        idx_wanted   = idx_ptp & idx_phase & ~idx_practice;
        
        % Some participants failed QC, so had no phase 2. Skip the loop for
        % them.
        if nnz(idx_wanted) == 0
            ctr = ctr + 1;
            continue
        end
        
        curr_data    = data_long(idx_wanted,:);
        
        % How many sessions were there?
        n_sessions = numel(unique(curr_data.session));
        
        for iSession = 1:n_sessions
            
            curr_session_data = curr_data(curr_data.session == iSession,:);
                       
            % Get the actual performance, correct/incorrect on each trial.
            % NaN means missed trials.
            y            = curr_session_data.correct;
            n_trials     = length(y);
            
            assert(n_trials == 42);
            
            %% Were the scores calculated right?
            
            % Reset variables
            rowidx_of_last_trial = [];
            n_trials_all_targets = [];
            
            % Which prompt toys were used? Get their names
%             [prompt_groups, prompt_names] = findgroups(curr_data.prompt_img_name);
            
            prompt_names = {'Sledge','Gingerman','Bear'};
            
            % Go trial by trial and check
            for iTrial = 1:n_trials
                iTrial
                curr_perf = curr_session_data.running_perf{iTrial};
                
                trunc_table = curr_session_data(1:iTrial,:);
                
                % For each toy check if the performance was right
                for iGroup = 1:numel(prompt_names)
                    
                    curr_toy = prompt_names{iGroup};
                    
                    % Get this toys performance
                    iPerf = trunc_table.correct(strcmp(...
                        trunc_table.prompt_img_name,curr_toy));
                    
                    % Change NaN's in iPerf to 0s
                    iPerf(isnan(iPerf)) = 0;
                    
                    % mean performance?
                    iMeanPerf = nanmean(iPerf);
                    
                    % Percentage performance?
                    iPercPerf{iPhase}{iSession}(iTrial,iGroup)  = iMeanPerf * 100;
                    
                    if isnan(iPercPerf{iPhase}{iSession}(iTrial,iGroup))
                        iPercPerf{iPhase}{iSession}(iTrial,iGroup) = 0;
                    end
                    
                    % Matches?
                    assert(iPercPerf{iPhase}{iSession}(iTrial,iGroup) == curr_perf(iGroup));
                end % iGroup                
            end % iTrial
        end % iSession
    end % iPhase
end % iPtp

%% Subfunctions
function n_trials_took = sliding_window(data,w,crit)
 
    % Pad the array in the beginning with 0.5
    data_padded = padarray(data,w/2,0.5,'pre');
    
    % Apply the window
    data_smoothed = movmean(data_padded,[w-1 0],'omitnan');

    data_smoothed_unpad = data_smoothed(w/2+1 : end);
    
    % How many trials did it take to reach criterion?
    n_trials_took = find(data_smoothed_unpad >= crit,1);
    
    if isempty(n_trials_took)
        n_trials_took = NaN;
    end
end