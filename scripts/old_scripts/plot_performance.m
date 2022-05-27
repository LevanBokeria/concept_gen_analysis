%% Draft script to analyze study 2 data

% This script plots the performances in separate sessions.
% Can filter whether you want the plot to be about all toys, or specific
% one.


clear; clc; close all;

dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

addpath(genpath(home));

task = 'pilots';

participants = table;

participants.oldStyle = [1,1,1,0,0,0,0,0]';

writeToExcel   = 0;
plotRThist     = 0;
savePerfFigure = 0;

for iPtp = 1:height(participants)
    
    % Reset the names variable
    names = {};
    
    curr_ptp = participants.names{iPtp};
    
    curr_style = participants.oldStyle(iPtp);
    
    %% Read and decode
    rawtext = fileread(fullfile(home,'data',task,['jatos_results_' curr_ptp '.txt']));

    % Split at the first component
    splitter_start = '[data_submission_start';
    [subject_split,matches] = split(rawtext,[splitter_start '---']);

    splitter_end = 'data_submission_end]';
    [split2,matches2] = split(subject_split{2},['---' splitter_end]);

    decoded = jsondecode(split2{1});

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
                
                % Change rt column to cell
                if isa(ses_tbl.rt(1),'double')
                    ses_tbl.rt = ...
                        num2cell(ses_tbl.rt);
                end
                
                if isa(ses_tbl.key_press(1),'double')
                    ses_tbl.key_press = ...
                        num2cell(ses_tbl.key_press);
                end               

                if isa(ses_tbl.correct(1),'double')
                    ses_tbl.correct = ...
                        num2cell(ses_tbl.correct);
                end                               
                
                phase_results = [phase_results; ses_tbl];

            end
       
        end
    else
        % Data coded new way, in the row format 
        phase_results = struct2table(decoded.outputData.phase_results);
        
        phase_results.ptp = cell(height(phase_results),1);
        [names{1:height(phase_results)}] = deal(curr_ptp);
        
        phase_results.ptp = names';
    end

    %% Record the participant table
    participants.data{iPtp} = phase_results;
    
end % for iPtp

%% Start plotting perf

% Should we filter by toy?
toy = 'Bear';

for iPtp = 1:height(participants)
    
    figure
    set(gcf,'Position',[    1.9557   -0.0677    1.1980    0.4200]*1000);
    
    curr_ptp = participants.names{iPtp}
    
    curr_data = participants.data{iPtp};
    
    % If the "correct" has empty cells, turn them to 0
    idx_empty = find(cellfun(@isempty,curr_data.correct));
    if ~isempty(idx_empty)
        curr_data.correct(idx_empty) = {0};
    end
    
    % how many phases?
    n_phases = numel(unique(curr_data.phase));
    
    for iPhase = 1:n_phases
        
        phase_table = curr_data(curr_data.phase == iPhase,:);
        
        if isa(phase_table.session(1),'cell')
            non_practice_idx = ~strcmp(phase_table.session,'practice');
            phase_table = phase_table(non_practice_idx,:);
            
            phase_table.session = cell2mat(phase_table.session);
        end
        
        n_sessions = numel(unique(phase_table.session));
       
        for iSes = 1:n_sessions
            
            session_table = phase_table(phase_table.session == iSes,:);
            subplot(3,6,(iPhase-1)*6 + iSes)
            
            if ~strcmp(toy,'all_toys')
                idx_curr_toy = strcmp(session_table.prompt_img_name,toy);
                session_table = session_table(idx_curr_toy,:);
            end
            
            switch toy
                case 'all_toys'
                    line_color = 'b';
                case 'Sledge'
                    line_color = 'r';
                case 'Gingerman'
                    line_color = [255,215,0]/255;
                case 'Bear'
                    line_color = 'g';
            end

            
            plot(cell2mat(session_table.correct),'Color',line_color,...
                'LineWidth',1.5);
            ylim([-0.1 1.1]);
            
            title(['Session ' int2str(iSes)]);
        end % iSes
        
    end    % iPhase
    
    sgtitle(['Ptp ' int2str(iPtp) ' All toys']);

    if savePerfFigure
        % Save the figure
        saveas(gcf,fullfile(home,'analysis',['ptp_' int2str(iPtp) '_' toy '.png']));
        
    end
end % iPtp

