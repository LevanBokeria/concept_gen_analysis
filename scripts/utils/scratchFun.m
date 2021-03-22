%% Dealing with incomplete participants

clear; clc;
dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

addpath(genpath(home));

% Suppress warnings
warning('off','MATLAB:table:RowsAddedExistingVars')

% Name of the file containing participant data from JATOS
jatos_filename_short = 'wave_3';
jatos_filename_full  = ['jatos_results_' jatos_filename_short '.txt'];

% Flags
saveCurrPtpVars   = 0; % Save the preprocessed files for the current participants?
plotPtpHist       = 0;
saveWaveData      = 0; % Save the combined feedback data for this wave, so we can manually look through it.

%% Read the filename and determine the number of participants

rawtext = fileread(fullfile(home,'data',jatos_filename_full));

% Split at the first component
[ptp_split,~] = split(rawtext,'[get_pid_comp_start---');

nPtp = numel(ptp_split) - 1;

for iPtp = 1:nPtp

    
    % Redefine ptp specific data
    qc_table = table;
    
    % Take this data
    curr_cell = ptp_split{iPtp+1};
    
    % Get the data submission section
    splitter_start = '[data_submission_start';
    [curr_cell,~] = split(curr_cell,[splitter_start '---']);
    
    % First, determine if this participant has a data submission component
    if numel(curr_cell) > 1
        assert(numel(curr_cell) == 2, 'too many data submissions!');
        ds_exists = 1;
    else
        ds_exists = 0;
    end
    
    if ~ds_exists
        
        
        
         % Separate out the consent component
        try
            
            splitter = '---';
            consent_split = split(curr_cell,splitter);
            splitter = '---]';
            [consent_split,matches] = split(consent_split{2},splitter);
            
            % Get the consent cell array only
            consent_struct = jsondecode(consent_split{1});
            
            % Get subject name
            curr_ptp  = consent_struct.prolific_ID;
            
        catch
            
            curr_ptp  = 'N/A'

        end 
        
        
        
        
    end


end
%%
% % %% Plot the learning script
% % clear; clc; 
% % 
% % close all;
% % 
% % dbstop if error;
% % 
% % trials = 0:50;
% % 
% % c = 0.1;
% % 
% % coeff = 0.5;
% % 
% % y_hat = 1 - coeff * exp(-c*trials);
% % 
% % 
% % plot(trials,y_hat);
% % grid on
% % xlim([-10 trials(end)]);
% % ylim([0 1.1]);
% % Draft script to analyze study 2 data
% clear; clc; close all;
% 
% dbstop if error;
% 
% % Check you're in the right directory 
% home = pwd;
%    
% [~,name,~] = fileparts(home);
% 
% if ~strcmp('concept_gen_analysis',name)
%     error('please change working directory to ./con_learn/concept_gen_analysis/');
% end
% 
% addpath(genpath(home));
% 
% task = 'pilots';
% 
% % jatosFileName = '5c6092ac4cc94200016a61be'; % failed at s3
% % jatosFileName = 'pilot25_int_fb';
% jatosFileName = 'pilot26';
% 
% % participants.names = {'5f0df9475f6bfa0f1288e014','5ece4950c4b2411f81c16674',...
% %                       '5e99f1ab02672c06c51487bd','5cc7650d4bb96d00018a95fd',...
% %                       '5f535ca509e9696557f043d6'};
% 
% oldStyle = 0;
% % participants.oldStyle     = [0,0,0,0,0];
% writeToExcel = 0;
% 
% %% Read and decode
% rawtext = fileread(fullfile(home,'data',task,['jatos_results_' jatosFileName '.txt']));
% 
% % Split at the first component
% splitter_start = '[data_submission_start';
% [subject_split,matches] = split(rawtext,[splitter_start '---']);
% 
% splitter_end = 'data_submission_end]';
% [split2,matches2] = split(subject_split{2},['---' splitter_end]);
% 
% % splitter_start = '[int_fb_start';
% % [subject_split,matches] = split(rawtext,[splitter_start '---']);
% % 
% % splitter_end = 'int_fb_end]';
% % [split2,matches2] = split(subject_split{2},['---' splitter_end]);
% 
% decoded = jsondecode(split2{1});
% 
% if oldStyle
%     % Data is coded old way, separate struct for each phase
%     
%     % Show me quick relevant data
%     qc_status = decoded.qc_status;
% 
%     % Running performance
%     running_perf = decoded.inputData.running_perf;
% 
%     %% Debriefs
%     int_debriefs   = table;
%     final_debriefs = table;
% 
%     for iPhase = 1:3
% 
%         phase_string = ['phase_' int2str(iPhase)];
% 
%         % Intermediate debriefs
%         if ~isempty(decoded.outputData.intermediate_feedback_results.(phase_string))
%             temp_int_debriefs_table = struct2table(decoded.outputData.intermediate_feedback_results.(phase_string));
%             phase_col = ones(height(temp_int_debriefs_table),1)*iPhase;
%             temp_int_debriefs_table.phase = phase_col;
% 
%             for iSes = 1:height(temp_int_debriefs_table)
% 
%                 temp_int_debriefs_table.session(iSes) = iSes;
%             end
% 
%             int_debriefs = [int_debriefs; temp_int_debriefs_table];
% 
%         end
% 
%         % Final debriefs
%         if ~isempty(decoded.outputData.([phase_string '_debrief_results']))
% 
%             temp_array = decoded.outputData.([phase_string '_debrief_results']){1,1};
% 
%             for iCell = 1:size(temp_array,1)
% 
%                 iTestPart = temp_array{iCell}.test_part;
% 
%                 if ismember(iTestPart,['general_questions_page','phase_1_debriefer',...
%                         'phase_2_debriefer','realize_mapping_page','phase_1_spatial_layout',...
%                         'phase_2_spatial_layout','phase_3_was_1_page'])
% 
%                     iResponse = jsondecode(temp_array{iCell}.responses);
% 
%                     final_debriefs.(iTestPart) = iResponse;
%                 end
%             end %iCell of final debriefer
%         end % if not isempty final debriefs
% 
%     end % iPhase
% 
% else
%     % Data coded new way, in the row format 
%     phase_results = struct2table(decoded.outputData.phase_results);
%     
%     % Intermediate debriefs:
%     int_debriefs = decoded.outputData.intermediate_feedback_results;
%     
%     % Final debriefs
%     if ~isempty(decoded.outputData.phase_1_debrief_results)
%         final_debriefs = decoded.outputData.phase_1_debrief_results;
%     elseif ~isempty(decoded.outputData.phase_2_debrief_results)
%         final_debriefs = decoded.outputData.phase_2_debrief_results;
%     elseif ~isempty(decoded.outputData.phase_3_debrief_results)
%         final_debriefs = decoded.outputData.phase_3_debrief_results;
%     end
%     
%     final_debriefs = final_debriefs{1}(2:end-1);
%     fb_tbl = table;
%     for iCell = 1:numel(final_debriefs)
%         ifb_decoded = struct2cell(jsondecode(final_debriefs{iCell}.responses));
%         
%         table_height = height(fb_tbl);
%         for iSubCell = 1:numel(ifb_decoded)
%             
%             fb_tbl.test_part{table_height+iSubCell} = final_debriefs{iCell}.test_part;
%             fb_tbl.responses{table_height+iSubCell} = ifb_decoded{iSubCell};
%         end
%     end
%     
%     %% Plot RT histograms
%     n_phases = numel(unique(phase_results.phase));
%     
%     for iPhase = 1:n_phases
%         
%         phase_table = phase_results(phase_results.phase == iPhase,:);
%         
%         % Remove practice trials
%         if isa(phase_table.session(1),'cell')
%             non_practice_idx = ~strcmp(phase_table.session,'practice');
%             phase_table = phase_table(non_practice_idx,:);
%             
%             phase_table.session = cell2mat(phase_table.session);
%         end
%         
%         % Find session == -1, those are also practice trials.
%         non_practice_idx = phase_table.session ~= -1;
%         phase_table = phase_table(non_practice_idx,:);
%         
%         n_sessions = numel(unique(phase_table.session));
%        
%         for iSes = 1:n_sessions
%             
%             session_table = phase_table(phase_table.session == iSes,:);
%             
%             rt = cell2mat(session_table.rt);
%             figure
%             histogram(rt,20);
%             ylim([0 15]);
%             xlim([0 11000]);
%             
%             title(['Subj ' jatosFileName ', Phase ' int2str(iPhase) ', session ' int2str(iSes)]);
% %             set(gcf,'Position',[21+((iSes-1)*300)   720-(iPhase-1)*330   300 216]);
%         end
%         
%     end % iPhase
%     %% Write to excel
%     if writeToExcel
%         % Change some of the columns, so it can be writte in the excel file
%         for iCell = 1:height(phase_results)
%             
%             phase_results.item_img_names{iCell} = phase_results.item_img_names{iCell}';
%             phase_results.item_img_paths{iCell} = phase_results.item_img_paths{iCell}';
%         end
%         
%         writetable(phase_results,fullfile(home,'analysis','pilots.xlsx'),...
%             'Sheet',jatosFileName)
%     end
%       
% end
