%% Draft script to analyze study 2 data

% This script will plot performances of participants continuously. Each
% phase gets one subplot. Each session start is demarcaded.
% Flags allow choosing to plot all participants together, or in separate
% figures. Also, whether to plot all toys data or specific toys.
% Additional parameters to choose how to smooth data

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
                  
participants.names = {'pilot20','pilot21','pilot22','pilot23','pilot24','pilot26','pilot27'}' 
% participants.names = {'pilot17','pilot18','pilot19'}' 
% participants.names = {'pilot11','pilot15','pilot16'}' 

% participants.oldStyle     = [1,1,zeros(1,6)]';
participants.oldStyle     = zeros(7,1);

nPtp = height(participants);

writeToExcel   = 0;
plotRThist     = 0;

plotOverlaid   = 1; % should participant lines be plotted together? Or separate figures for each ptp?
savePerfFigure = 0;

averageAcrossPtp = 0;

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

    participants.congruency(iPtp) = decoded.inputData.congruency;
    
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

% Figure position second screen
% fig_pos = [1.9130   -0.1257    1.2800    0.6073]*1000; 
% fig_pos = [0.0010    0.0410    1.2800    0.6073]*1000;
% Single screen
fig_pos = [0.0250    0.1470    1.2147    0.4200]*1000; 

% Should we filter by toy?
% toys = {'Sledge','Gingerman','Bear'};
toys = {'all_toys'};

% Smoothing parameters
smooth_methods = {'none','gaussian'};
% smooth_methods = {'none'};
window_sizes = [5];


if averageAcrossPtp
    
    saveFolder = fullfile('smoothing','all_toys');
    saveLoc    = fullfile(home,'analysis',saveFolder);
    if ~exist(saveLoc)
        mkdir(saveLoc);
    end
    
    figure
    set(gcf,'Position',fig_pos)
    
    % Create a cell array to hold participant data, such that we can later
    % average across participants, despite different lenghts of trials
    avg_ptp_data = {};
    
    for iPhase = 1:3
        
        subplot(3,1,iPhase);
        title(['Phase ' int2str(iPhase)]);

        % Plot demarcation lines
        for iLine = 1:6
            x_coord = 1 + (iLine-1)*42;
            line([x_coord,x_coord],[1.2,1.5],'Color','k',...
                'LineWidth',1,'LineStyle','-');
            hold on
        end
        
        % Calculate the average across participants
        for iPtp = 1:nPtp
        
            curr_ptp = participants.names{iPtp}
            
            curr_data = participants.data{iPtp};
            
            % If the "correct" has empty cells, turn them to 0
            if isa(curr_data.correct(1),'cell')
                idx_empty = find(cellfun(@isempty,curr_data.correct));
                if ~isempty(idx_empty)
                    curr_data.correct(idx_empty) = {0};
                end
                curr_data.correct = cell2mat(curr_data.correct);
            end
            
            phase_table = curr_data(curr_data.phase == iPhase,:);
            
            % Remove practice trials
            if isa(phase_table.session(1),'cell')
                non_practice_idx = ~strcmp(phase_table.session,'practice');
                phase_table = phase_table(non_practice_idx,:);
                
                phase_table.session = cell2mat(phase_table.session);
            end
            
            % Find session == -1, those are also practice trials.
            non_practice_idx = phase_table.session ~= -1;
            phase_table = phase_table(non_practice_idx,:);
            
            avg_ptp_data{iPtp,iPhase} = phase_table.correct;
        end % iPtp
        
        % Find the one with most number of trials
        sizes = cellfun(@size,avg_ptp_data(:,iPhase),'UniformOutput',false);
        sizes = cell2mat(sizes);
        
        [max_val,idx] = max(sizes(:,1));
        
        iPhaseNan{iPhase} = nan(nPtp,max_val);
        
        % Populate the nan table with participant values
        for iiPtp = 1:nPtp
            
            % How many values are there now?
            nVals = length(avg_ptp_data{iiPtp,iPhase});
            
            iPhaseNan{iPhase}(iiPtp,1:nVals) = avg_ptp_data{iiPtp,iPhase};
            
        end
        
        % Take the average
        iPhaseNanMean{iPhase} = nanmean(iPhaseNan{iPhase},1);

        % Plot this average
         plot(iPhaseNanMean{iPhase},...
             'LineWidth',1);
         ylim([-0.1 1.5]);
         xlim([0 42*6]);
         xticks(0:5:42*6);
         a = get(gca,'XTickLabel');
         set(gca,'XTickLabel',a,'FontName','Times','fontsize',7)
        
    end % iPhase 

    sgtitle('Averaged participants');
    
    saveas(gcf,fullfile(saveLoc,'avg_ptps_all_toys.png'));
end    
% else
    % So we are not averaging over participants
    
    for iToy = 1:numel(toys)

        curr_toy = toys{iToy};

        % Set the color of the line
        switch curr_toy
            case 'all_toys'
                line_color = 'b';
            case 'Sledge'
                line_color = 'r';
            case 'Gingerman'
                line_color = [255,215,0]/255;
            case 'Bear'
                line_color = 'g';
        end

        saveFolder = fullfile('smoothing','prolific',curr_toy);
        saveLoc    = fullfile(home,'analysis',saveFolder);
        if ~exist(saveLoc)
            mkdir(saveLoc);
        end

        for iMethod = 1:size(smooth_methods,2)

            curr_method = smooth_methods{iMethod};

            for iWS = 1:numel(window_sizes)

                curr_window_size = window_sizes(iWS)

    %             close all

                if plotOverlaid
                    figure
                    set(gcf,'Position',fig_pos)
                end

                for iPtp = 1:height(participants)

                    y_offset = iPtp*0.02;

                    if ~plotOverlaid
                        figure
                        set(gcf,'Position',fig_pos)
                        y_offset = iPtp*0.00;
                    end

                    curr_ptp = participants.names{iPtp}

                    curr_data = participants.data{iPtp};

                    % If the "correct" has empty cells, turn them to 0
                    if isa(curr_data.correct(1),'cell')
                        idx_empty = find(cellfun(@isempty,curr_data.correct));
                        if ~isempty(idx_empty)
                            curr_data.correct(idx_empty) = {0};
                        end
                        curr_data.correct = cell2mat(curr_data.correct);
                    end
                    % how many phases?
                    n_phases = numel(unique(curr_data.phase));

                    for iPhase = 1:2 %n_phases

                        phase_table = curr_data(curr_data.phase == iPhase,:);

                        % Remove practice trials
                        if isa(phase_table.session(1),'cell')
                            non_practice_idx = ~strcmp(phase_table.session,'practice');
                            phase_table = phase_table(non_practice_idx,:);

                            phase_table.session = cell2mat(phase_table.session);
                        end

                        % Find session == -1, those are also practice trials.
                        non_practice_idx = phase_table.session ~= -1;
                        phase_table = phase_table(non_practice_idx,:);

                        n_sessions = numel(unique(phase_table.session));                

                        subplot(2,1,iPhase)

                        if strcmp(curr_toy,'all_toys')

                            % So if its all toys, just plot the whole data

                            % First, plot the session demarcation lines
                            if iPtp == 1 | plotOverlaid == 0
                                for iLine = 1:4
                                    x_coord = 1 + (iLine-1)*42;
                                    line([x_coord,x_coord],[1.2,1.5],'Color','k',...
                                        'LineWidth',1,'LineStyle','-');
                                    hold on
                                end
                                %             hold off
                            end % if iPtp == 1

                            if isa(phase_table.correct(1),'cell')
                                data = cell2mat(phase_table.correct)';
                            else
                                data = phase_table.correct';
                            end
                            data_smoothed = smooth_performance(data,curr_window_size,curr_method);

                            plot(data_smoothed + y_offset,...
                                'LineWidth',1);
                            ylim([-0.1 1.5]);
                            xlim([0 42*4]);
                            xticks(0:5:42*4);
                            a = get(gca,'XTickLabel');
                            set(gca,'XTickLabel',a,'FontName','Times','fontsize',7)
                        else

                            % So we're plotting just the data for specific toys.

                            idx_curr_toy = strcmp(phase_table.prompt_img_name,curr_toy);
                            phase_table = phase_table(idx_curr_toy,:);

                            x_coord = 1;

                            data = [];

                            nTrialsPerSession = 14;

                            for iSes = 1:n_sessions

                                session_table = phase_table(phase_table.session == iSes,:);

                                % Get the data and truncate it to 13 trials 
                                data = [data session_table.correct(1:nTrialsPerSession)'];

                                x_coord = 1 + (iSes-1)*nTrialsPerSession;

                                line([x_coord,x_coord],[1.2,1.5],'Color','k',...
                                    'LineWidth',1,'LineStyle','-');
                                hold on
                            end

                            data_smoothed = smooth_performance(data,curr_window_size,curr_method);

                            if plotOverlaid
                                plot(data_smoothed + y_offset,...
                                    'LineWidth',1);                        
                            else
                                plot(data_smoothed + y_offset,...
                                    'LineWidth',1,'Color',line_color);
                            end
                            ylim([-0.15 1.5]);
                            xlim([0 nTrialsPerSession*6]); % so if a toy has 14 trials, times 6 sessions

                            xticks(0:nTrialsPerSession*6);
                            a = get(gca,'XTickLabel');
                            set(gca,'XTickLabel',a,'FontName','Times','fontsize',7)
                        end % if all_toys or other toy

%                         title(['Phase ' int2str(iPhase)]);

                        if plotOverlaid
                            hold on
                        end        

                        grid on
                    end    % iPhase

                    if ~plotOverlaid

                        sgtitle(['ptp' int2str(iPtp) '_' curr_toy '. Smooth method "' ...
                            curr_method '", window=' int2str(curr_window_size), ...
                            ', congruent = ' int2str(participants.congruency(iPtp))],...
                            'Interpreter','none');

                        if savePerfFigure
                            % Save the figure
                            saveas(gcf,fullfile(saveLoc,['ptp_' int2str(iPtp) '_' ...
                                curr_toy '_' curr_method '_wind_' int2str(curr_window_size) '.png']));

                        end
                    end
                end % iPtp

                if plotOverlaid
%                     sgtitle(['All ptps, ' curr_toy '. Smooth method "' curr_method ...
%                         '", window=' int2str(curr_window_size)],...
%                         'Interpreter','none');

                    if savePerfFigure

                        saveas(gcf,fullfile(saveLoc,['all_ptp_' curr_toy '_' ...
                            curr_method '_wind_' int2str(curr_window_size) '.png']));

                    end
                end
            end % for iWindowSize
        end % for iMethod

    end % for iToy
% end
function [data_out] = smooth_performance(data_in,window_length,smooth_method)
    
    if (nargin < 2)
        window_length = 5;
        smooth_method = 'movmean';
    end
    
    %% Pad the data at the ends
    pad_value = 0.5;
    pad_amount = (window_length-1)/2;
    
    data_in_padded = padarray(data_in,[0,pad_amount],pad_value,'pre');
    
    %% Smooth the data
    if (strcmp(smooth_method,'none'))
        data_in_padded_smoothed = data_in_padded;
    else
        data_in_padded_smoothed = smoothdata(data_in_padded,smooth_method,window_length);
    end
    %% Remove padding
    data_in_smoothed_depadded = data_in_padded_smoothed(pad_amount+1:end);

    %% Output data
    data_out = data_in_smoothed_depadded;

end