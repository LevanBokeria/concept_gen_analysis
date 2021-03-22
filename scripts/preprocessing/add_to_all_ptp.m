function add_to_all_ptp(batch_number,wave_name,overwriteData,saveCombinedData)

%% If not arguments, load wave_1 data
if (nargin < 1)
    batch_number    = 'batch_0';
    wave_name       = 'wave_5';
    overwriteData    = 1;
    saveCombinedData = 1;
end

% Check you're in the right directory 
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

batch_wave_name = [batch_number '_' wave_name];

%% Load the data for all ptps:
[long_form_data_all_ptp,results_table_all_ptp,...
    debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
    load_all_participant_files(home);

%% Load the data for this wave
wave_data = load(fullfile(home,'results','analysis','waveData',...
    batch_wave_name,'wave_results_debrief_and_feedback.mat'));
wave_long_form_data = load(fullfile(home,'results','analysis',...
    'waveData',batch_wave_name,'wave_long_form_data.mat'));
wave_long_form_data = wave_long_form_data.wave_long_form_data;


wave_results      = wave_data.wave_results;
wave_debriefing   = wave_data.wave_debriefing;
wave_int_feedback = wave_data.wave_int_feedback;

%% Check if any participant from this wave is already in the data
for iPtp = 1:height(wave_results)
    curr_ptp = wave_results.ptp{iPtp};
    
    % For long form data
    long_form_data_all_ptp = ...
        process_participants(long_form_data_all_ptp,wave_long_form_data,curr_ptp);
    
    % For results table
    results_table_all_ptp = ...
        process_participants(results_table_all_ptp,wave_results,curr_ptp);

    % For debriefs
    debrief_tbl_all_ptp = ...
        process_participants(debrief_tbl_all_ptp,wave_debriefing,curr_ptp);
    
    % For fb intermediate table
    fb_inter_tbl_all_ptp = ...
        process_participants(fb_inter_tbl_all_ptp,wave_int_feedback,curr_ptp);
    
end

%% Save the combined data?
if saveCombinedData
    
    fprintf('Saving comibned data...\n');
    
    save(fullfile(home,'results','analysis',...
        'long_form_data_all_ptp.mat'),'long_form_data_all_ptp');
    save(fullfile(home,'results','analysis',...
        'results_table_all_ptp.mat'),'results_table_all_ptp');
    save(fullfile(home,'results','analysis',...
        'debrief_tbl_all_ptp.mat'),'debrief_tbl_all_ptp');
    save(fullfile(home,'results','analysis',...
        'fb_inter_tbl_all_ptp.mat'),'fb_inter_tbl_all_ptp');
    
    fprintf('Saved all the data!\n');
end

%% Subfunction
    function all_ptp_data = process_participants(all_ptp_data,wave_data,curr_ptp)
        
        %% For long form data
        
        % Get this participants data
        idx = strcmp(wave_data.ptp,curr_ptp);
        curr_ptp_data = wave_data(idx,:);
        
        % Is this participants data not in the all ptp data?
        if height(all_ptp_data) == 0 || ...
                nnz(strcmp(all_ptp_data.ptp,curr_ptp)) == 0
            
            % So data is not in the long form data, add it
            if height(all_ptp_data) == 0
                all_ptp_data = curr_ptp_data;
            else
                all_ptp_data = [all_ptp_data; curr_ptp_data];
            end
            
        elseif overwriteData
            
            % So it is in the all_ptp data, but we want to overwrite
            
            fprintf(['Ptp ' curr_ptp ' already in the all ptp long_form data\n']);
            fprintf('Overwriting...\n\n');
            
            % Find indices of this ptp in all_ptp data
            idx_in_all_data = strcmp(all_ptp_data.ptp,curr_ptp);
            
            % Sanity check
            assert(nnz(idx_in_all_data) == height(curr_ptp_data));
            
            % Change this rows with data from the wave
            all_ptp_data(idx_in_all_data,:) = curr_ptp_data;
            
        else
            % So its in the all_ptp, and don't overwrite
            fprintf(['Ptp ' curr_ptp ' already in the all ptp long_form data\n']);
            fprintf('NOT overwriting...\n\n');
        end
        
        
    end

end