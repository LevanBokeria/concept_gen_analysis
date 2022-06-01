function preprocess_wave_data_submission_module()

%% Description

% Process just the data submission modules


%% Define global variables and conditions
clear; close all;

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
batch_number    = 'batch_3';
wave_name       = 'wave_17_ds';
batch_wave_name = [batch_number '_' wave_name];

% Flags
saveCurrPtpVars   = 1; % Save the preprocessed files for the current participants?
saveWaveData      = 1; % Save the combined feedback data for this wave, so we can manually look through it.

%% Read the filename and determine the number of participants

rawtext = fileread(fullfile(home,...
    'data','jatos_data_submission_exports',...
    ['jatos_results_' batch_wave_name '.txt']));

% Split at the first component
[ptp_split,~] = split(rawtext,'[data_submission_start---');

nPtp = numel(ptp_split) - 1;

%% Start the loop to preprocess

% Predefine the tables for this wave
wave_results        = table;
wave_debriefing     = table;
wave_int_feedback   = table;
wave_long_form_data = table;

% Start the ptp loop
for iPtp = 1:nPtp
    iPtp
    % Reset variables
    curr_ptp              = [];
    idebrief_wide_format  = table;
    ifb_int_long_format   = table;
    results_table         = table;
    long_form_data        = table;
    no_ds_with_debriefing = 0;
    
    % Take this data
    curr_cell = ptp_split{iPtp+1};
    
    % Get the data submission section
%     splitter_start = '[data_submission_start';
%     [curr_cell,~] = split(curr_cell,[splitter_start '---']);
%     
%     % First, determine if this participant has a data submission component
%     if numel(curr_cell) > 1
%         assert(numel(curr_cell) == 2, 'too many data submissions!');
%         ds_exists = 1;
%     else
%         ds_exists = 0;
%     end
    ds_exists = 1;
    
    if ds_exists
        
        splitter_end = 'data_submission_end]';
        [split2,~] = split(curr_cell,['---' splitter_end]);
        
        decoded = jsondecode(split2{1});

        %% Get the variables to record
        curr_ptp                = clean_name(decoded.prolific_ID);
        iCongruency             = decoded.inputData.congruency;
        iQCTable                = struct2table(decoded.qc_status);
        iConcept_phase_1        = decoded.inputData.concepts.phase_1.concept_space;
        iConcept_phase_2        = decoded.inputData.concepts.phase_2.concept_space;
        iArr_phase_1            = decoded.inputData.basic_parameters.targetPoints.phase_1;
        iArr_phase_2            = decoded.inputData.basic_parameters.targetPoints.phase_2;
        
        targetNamesUsed_phase_1 = decoded.inputData.basic_parameters.targetNamesUsed.phase_1;
        targetNamesUsed_phase_2 = decoded.inputData.basic_parameters.targetNamesUsed.phase_2;
        
        iProgressState          = decoded.progress_state;
        
        % At what stage did the experiment end?
        tempPhaseResults = struct2table(decoded.outputData.phase_results);
        iMaxPhase   = tempPhaseResults.phase(end);
        iMaxSession = tempPhaseResults.session(end);
        
%         % How long did they take to read instructions?
%         iPhase_1_instr_table = struct2table(jsondecode(decoded.outputData.instructions_results.phase_1.view_history));
%         iPhase_1_instr_dur   = sum(iPhase_1_instr_table.viewing_time);
    else
        
        splitter = '---';
        components_split = split(curr_cell,splitter);
        
%         input(['This ptp had no data_submission, trying to decode...']);
        
        nCells = numel(components_split);
        decoded = jsondecode(components_split{nCells-1});
        
        curr_ptp                = clean_name(decoded.prolific_ID)
        
        try
            iCongruency             = decoded.inputData.congruency;
            iQCTable                = struct2table(decoded.qc_status);
            iConcept_phase_1        = decoded.inputData.concepts.phase_1.concept_space;
            iConcept_phase_2        = decoded.inputData.concepts.phase_2.concept_space;
            iArr_phase_1            = decoded.inputData.basic_parameters.targetPoints.phase_1;
            iArr_phase_2            = decoded.inputData.basic_parameters.targetPoints.phase_2;
            
            targetNamesUsed_phase_1 = decoded.inputData.basic_parameters.targetNamesUsed.phase_1;
            targetNamesUsed_phase_2 = decoded.inputData.basic_parameters.targetNamesUsed.phase_2;
            
            iProgressState          = decoded.progress_state;
            
            % At what stage did the experiment end?
            if ~isempty(decoded.outputData.phase_results)
                
                % So we have at least some data
                tempPhaseResults = struct2table(decoded.outputData.phase_results);
                iMaxPhase   = tempPhaseResults.phase(end);
                iMaxSession = tempPhaseResults.session(end);
                
            else
                % If its empty, then didnt even complete the practice
                iMaxPhase   = NaN;
                iMaxSession = NaN;
                
            end
            
            % Maybe no data submission due to a bug, but still did the debriefing?
            no_ds_with_debriefing = no_ds_ptp_with_debriefing(curr_ptp,decoded);
            
        catch
            % So, ptp quit before conditions were even assigned
            iCongruency             = NaN;
            
            % Create an empty table with field names for iQCTable
            iQCTable                = cell2table(num2cell(NaN(1,8)), ...
                'VariableNames', {'global_pass', 'rt_pass', ...
                'uniform_resp_perc_pass', 'min_perf_pass', ...
                'practice_pass', 'perc_max_misses_pass', ...
                'max_training_sess_pass', 'min_time_instruct_pass'});
            
            iConcept_phase_1        = NaN;
            iConcept_phase_2        = NaN;
            iArr_phase_1            = NaN;
            iArr_phase_2            = NaN;
            
            targetNamesUsed_phase_1 = NaN;
            targetNamesUsed_phase_2 = NaN;
            
            iProgressState          = 'condition not assigned';
            
            iMaxPhase = NaN;
            iMasSession = NaN;
        end
    end
    
    
        %% Record metadata in the QC table
        
        % Put all of this in a table
        results_table.ptp{1}                     = curr_ptp;
        results_table.progress_state{1}          = iProgressState;
        results_table.maxPhase(1)                = iMaxPhase;
        results_table.maxSesssion(1)             = iMaxSession;

        % Record that the data submission component exists for this ptp
        results_table.data_submitted             = ds_exists;        
        
        % Add the full QC table to it
        results_table = horzcat(results_table,iQCTable);
        
        % Place to later store info on whether RT was above 3SD of group average
        results_table.phase_1_rt_qc_pass         = NaN;
        results_table.phase_2_rt_qc_pass         = NaN;
        
        % Place to later store info on whether debriefing or intermediate
        % feedback made them fail QC
        results_table.debrief_qc_pass            = NaN;
        results_table.fb_int_qc_pass             = NaN;
        
        % Get the details of the condition
        results_table.congruency                 = iCongruency;
        results_table.concept_phase_1{1}         = iConcept_phase_1;
        results_table.concept_phase_2{1}         = iConcept_phase_2;
        results_table.arr_phase_1{1}             = iArr_phase_1;
        results_table.arr_phase_2{1}             = iArr_phase_2;
        results_table.targetNamesUsed_phase_1{1} = targetNamesUsed_phase_1;
        results_table.targetNamesUsed_phase_2{1} = targetNamesUsed_phase_2;
        
        % Add placeholders for the intercept and exponent of the learning
        % rate
        results_table.learning_rate_intercept(1) = NaN;
        results_table.learning_rate_exponent(1)  = NaN;

        %% Clean any irregularities
        results_table = clean_irregularities(curr_ptp,results_table);
        
        
        %% Get the long form data across all sessions and all phases
        long_form_data = get_long_form_data(decoded,results_table);
        
        %% Calculate session 1, session 2, and averaged ses1:2 for each phase
        avg_perf = [];
        for iPhase = 1:2
            for iSes = 1:2
                
                try
                    
                    idx_phase     = long_form_data.phase == iPhase;
                    idx_ses       = long_form_data.session == iSes;
                    idx_phase_ses = idx_phase & idx_ses;
                    
                    avg_perf(iPhase,iSes) = nanmean(long_form_data.correct(idx_phase_ses));
                    
                catch
                   
                    % So, long_form_data doesn't even exist, probably
                    % because the ptp quit during instructions or practice
                    avg_perf(iPhase,iSes) = NaN;
                    
                end
                
                results_table.(['phase_' int2str(iPhase) '_ses_' int2str(iSes) ...
                    '_perf'])(1) = avg_perf(iPhase,iSes);
            end
            
            % Get average over the 2 sessions
            results_table.(['phase_' int2str(iPhase) '_ses_1_2_perf'])(1) = ...
                mean(avg_perf(iPhase,:));
        end
        
        %% Get intermediate debriefs
        
        % Intermediate debriefs:
        try
            ifb_int_long_format = struct2table(decoded.outputData.intermediate_feedback_results);
            
            ptp = {};
            [ptp{1:height(ifb_int_long_format)}] = deal(curr_ptp);
            ptp = ptp';
            
            ifb_int_long_format = addvars(ifb_int_long_format,ptp,'Before','rt');
           
        catch
            % If fails, so haven't even gotten to any intermediate debriefs
            % Must have quit before practice or right after
            i_int_fb_exists = 0;
        end

        
        %% Final debriefs
        idebrief_wide_format = get_final_debriefs(decoded);
                
        %% Save participant specific things
        if saveCurrPtpVars
            
            fprintf('Saving current ptp data...\n');
            
            if ~exist(fullfile(home,'results','preprocessed_data',curr_ptp))
                mkdir(fullfile(home,'results','preprocessed_data',curr_ptp));
            end
            
            % Save the decoded variable in the folder
            save(fullfile(home,'results','preprocessed_data',curr_ptp,...
                'decoded.mat'),'decoded')
            save(fullfile(home,'results','preprocessed_data',curr_ptp,...
                'long_form_data.mat'),'long_form_data')
            
            % Save the feedback as a mat tile
            save(fullfile(home,'results','preprocessed_data',curr_ptp,...
                'idebrief_wide_format.mat'),'idebrief_wide_format');
            
            % Save the intermediate debriefs
            save(fullfile(home,'results','preprocessed_data',curr_ptp,...
                'ifb_int_long_format.mat'),'ifb_int_long_format');
            
        end
        
        %% Add this to the wave specific data
        combine_wave_data(long_form_data,results_table,idebrief_wide_format,...
            ifb_int_long_format);
    
end % iPtp


%% Save wave data
if saveWaveData
    
    fprintf('Writing wave-specific data and excel files...\n');
    
    % Save these as mat files
    if ~exist(fullfile(home,'results','analysis','waveData',...
        batch_wave_name))
        mkdir(fullfile(home,'results','analysis','waveData',...
        batch_wave_name));
    end
    save(fullfile(home,'results','analysis','waveData',...
        batch_wave_name,'wave_results_debrief_and_feedback.mat'),...
        'wave_debriefing','wave_int_feedback','wave_results');
    save(fullfile(home,'results','analysis','waveData',...
        batch_wave_name,'wave_long_form_data.mat'),...
        'wave_long_form_data');
    
    % Save these as excel files
    writetable(wave_debriefing,...
        fullfile(home,'results','analysis','waveData',...
        batch_wave_name,'wave_debriefing.xlsx'));
    writetable(wave_int_feedback,...
        fullfile(home,'results','analysis','waveData',...
        batch_wave_name,'wave_int_feedback.xlsx'));    

    fprintf('Finished!\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subfunctions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function long_form_data = get_long_form_data(decoded,results_table)
        
        % Data coded new way, in the row format
        try
            long_form_data = struct2table(decoded.outputData.phase_results);
            % If missed, put NA there, for each column of data
            cols_to_touch = {'rt','key_press','correct'};
            
            for iCol = 1:numel(cols_to_touch)
                curr_col = cols_to_touch{iCol};
                
                if isa(long_form_data.(curr_col)(1),'cell') % so we have cells, theres a miss
                    idx_empty = find(cellfun(@isempty,long_form_data.(curr_col)));
                    if ~isempty(idx_empty)
                        long_form_data.(curr_col)(idx_empty) = {NaN};
                    end
                    
                    % First, turn any strings to numbers
                    if isa(long_form_data.(curr_col){1},'char')
                        long_form_data.(curr_col) = cellfun(@str2num,long_form_data.(curr_col),'UniformOutput',0);
                    end
                    % Now, turn the cell array to mat
                    long_form_data.(curr_col) = cell2mat(long_form_data.(curr_col));
                end
            end
            
            % Add the participant name to the table
            tbl_to_add = repmat(results_table,height(long_form_data),1);
            
            long_form_data = [tbl_to_add, long_form_data];
        catch
            
            % So no long form data. Must have quit before practice trials
            no_long_form_data = 1;
            
            long_form_data = table;
        end
    end

    function idebrief_wide_format = get_final_debriefs(decoded)
        
        idebrief_wide_format = table;
        
        % Add the ptp name to it
        idebrief_wide_format.ptp{1} = curr_ptp;
        idebrief_wide_format.data_submitted(1) = ds_exists;
        try
            idebrief_wide_format.global_pass(1)    = decoded.qc_status.global_pass;
        catch
            try
                idebrief_wide_format.global_pass(1)    = decoded.qc_status.global_pass;
            catch
                idebrief_wide_format.global_pass(1)    = NaN;
            end
        end
        
        if ds_exists | no_ds_with_debriefing
            
            % Reset end_point
            debrief_end_point = [];
            
            % Final debriefs
            if ~isempty(decoded.outputData.phase_1_debrief_results)
                iFinal_debriefs = decoded.outputData.phase_1_debrief_results;
                % Record where the experiment ended
                debrief_end_point = 'phase_1';
            elseif ~isempty(decoded.outputData.phase_2_debrief_results)
                iFinal_debriefs = decoded.outputData.phase_2_debrief_results;
                debrief_end_point = 'phase_2';
            elseif ~isempty(decoded.outputData.phase_3_debrief_results)
                iFinal_debriefs = decoded.outputData.phase_3_debrief_results;
                debrief_end_point = 'phase_3';
            end
            
            % Put the debriefs in a nice table
            
            if isa(iFinal_debriefs,'cell')
                
                % So practice stage has been passed.
                iFinal_debriefs = iFinal_debriefs{1}(2:end-1);
            for iCell = 1:numel(iFinal_debriefs)
                ifb_decoded = struct2cell(jsondecode(iFinal_debriefs{iCell}.responses));
                ifb_page_name = iFinal_debriefs{iCell}.test_part;
                
                for iSubCell = 1:numel(ifb_decoded)
                    if isempty(ifb_decoded{iSubCell})
                        ifb_decoded{iSubCell} = NaN;
                    end
                    idebrief_wide_format.([ifb_page_name '_Q' ...
                        int2str(iSubCell-1)]){1} = ...
                        ifb_decoded{iSubCell};
                end
            end
            
            else
                % These are debriefs obtained after practice qc failure.
                % They are in the format of a struct, not a cell array of
                % structs
                debrief_end_point = 'practice_trials';
            end
            
            % Where did this participant end the experiment?
            idebrief_wide_format.debrief_end_point{1} = debrief_end_point;            
            
        end % if ds_exists
    end % idebrief_wide_format

    function combine_wave_data(long_form_data,results_table,idebrief_wide_format,...
            ifb_int_long_format)

        % Check if this is the first participant
        if height(wave_long_form_data) == 0
            wave_long_form_data = long_form_data;
        else
            wave_long_form_data = [wave_long_form_data; long_form_data];
        end
        
        % Now for QC table
        if height(wave_results) == 0
            wave_results = results_table;
        else
            wave_results = [wave_results; results_table];
        end
        
        % Now for the final debriefs
        
        if height(wave_debriefing) == 0
            % Make a fresh one, with field names predefined.
            wave_debriefing                              = table;
            wave_debriefing.ptp                      {1} = {};
            wave_debriefing.data_submitted           (1) = NaN;
            wave_debriefing.global_pass              (1) = NaN;
            wave_debriefing.debrief_end_point        {1} = {};
            wave_debriefing.debrief_qc_pass          (1) = NaN;
            wave_debriefing.general_questions_page_Q0{1} = {};
            wave_debriefing.general_questions_page_Q1{1} = {};
            wave_debriefing.general_questions_page_Q2{1} = {};
            wave_debriefing.general_questions_page_Q3{1} = {};
            wave_debriefing.general_questions_page_Q4{1} = {};
            wave_debriefing.phase_1_debriefer_Q0     {1} = {};
            wave_debriefing.phase_1_debriefer_Q1     {1} = {};
            wave_debriefing.phase_2_debriefer_Q0     {1} = {};
            wave_debriefing.phase_2_debriefer_Q1     {1} = {};
            wave_debriefing.realize_mapping_page_Q0  {1} = {};
            wave_debriefing.phase_1_spatial_layout_Q0{1} = {};
            wave_debriefing.phase_2_spatial_layout_Q0{1} = {};
        end 
              
        % Add to the table
        all_fields = wave_debriefing.Properties.VariableNames;
        
        % If the table is fresh, start at the first row, else start after
        % the last:
        if isempty(wave_debriefing.ptp{1})
            startRow = 1;
        else
            startRow = height(wave_debriefing)+1;
        end
        
        for iF = 1:numel(all_fields)
            
            curr_field = all_fields{iF};
            
            % If this field exists for this participant, add the data
            if nnz(strcmp(idebrief_wide_format.Properties.VariableNames,...
                    curr_field)) == 1
                
                if isa(wave_debriefing.(curr_field)(1),'cell')
                    wave_debriefing.(curr_field){startRow} = ...
                        idebrief_wide_format.(curr_field){1};
                else
                    wave_debriefing.(curr_field)(startRow) = ...
                        idebrief_wide_format.(curr_field)(1);
                end
                
            else
                % So it doesn't exist for this participant. They must
                % have ended early. So add NaN
                try
                    wave_debriefing.(curr_field){startRow} = NaN;
                catch
                    wave_debriefing.(curr_field)(startRow) = NaN;
                end
            end
        end
        
        % Now for intermediate debriefs
        if height(wave_int_feedback) == 0
            wave_int_feedback = ifb_int_long_format;
        else
            wave_int_feedback = [wave_int_feedback; ifb_int_long_format];
        end
    end % combine wave data

    function name_out = clean_name(name_in)
       
        name_worked  = name_in;
        empty_exists = 1;
        counter      = 0;
        
        while empty_exists
            if name_worked(end) == ' '
                fprintf(['Empty char in ID for ' name_in '\n']);
                name_worked = name_worked(1:end-1);
            else
                empty_exists = 0;
            end
            counter = counter + 1;
        end % while loop
        name_out = name_worked;
    end

    function [data_out] = clean_irregularities(curr_ptp,data_in)
        
        % For clarify, assign a separate variable for worked on data.
        data_worked_on = data_in;
        
        if strcmp(curr_ptp,'5f561a95aa1c4ea13672f138') | ...
                strcmp(curr_ptp,'5fcfa3168335430d143b4431')
            
            % This participant missed too many trials. In the JS, I had a
            % typo where I created a new object called
            % "perc_max_missed_pass", instead of assigning a 0 to an
            % existing object called "perc_max_misses_pass". This the
            % qc_status object had 1 extra field, and thus the
            % results_table has 1 extra field. Fix this issue.
            data_worked_on.perc_max_misses_pass = data_worked_on.perc_max_missed_pass;
            data_worked_on.perc_max_missed_pass = [];
        end
        
        data_out = data_worked_on;
                
    end

    function [no_ds_with_debriefing] = no_ds_ptp_with_debriefing(curr_ptp,decoded)
        
        % This function will check if the current participant has no data
        % submission module, but still has debriefing completed. This would
        % happen sometimes if for some reason, after the debriefing was
        % done, transferring the data to JATOS/prolific happened
        % incompletely. Typically, these ptps have "NOCODE" in their prolific
        % column.
        
        if ~isempty(decoded.outputData.phase_1_debrief_results) | ...
                ~isempty(decoded.outputData.phase_2_debrief_results)
            
            no_ds_with_debriefing = 1;
        else
            no_ds_with_debriefing = 0;
        end
    end

end % preprocess function