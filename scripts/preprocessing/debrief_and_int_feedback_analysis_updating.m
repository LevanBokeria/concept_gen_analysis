%% Description:
% 
% Once you've downloaded and preprocessed the latest wave data, and then 
% looked through their debriefing and intermediate feedback, and scored if 
% they passed or failed those QCs, run this script with the updated 
% variables still in the working memory.
% 
% This script will read the wave_results table and see which ptps passed qc for 
% debriefing and intermediate feedback. It will then copy that information 
% to the wave_debriefing file itself.
clear; clc;
dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

%% Read the files
batch_number    = 'batch_3';
wave_name       = 'wave_16_ds';
batch_wave_name = [batch_number '_' wave_name];

load(fullfile(home,'results','analysis','waveData',...
    batch_wave_name,'wave_results_debrief_and_feedback.mat'));
load(fullfile(home,'results','analysis','waveData',...
    batch_wave_name,'wave_long_form_data.mat'));

%% Start updating 

% Have you scored the debriefing?
input('Scored the debriefing in the excel file?\n')

% Have you copied the scoring to the mat file in working memory?
input('Have you copied the scoring to the wave_results mat file in working memory?\n');

% Have you scored the intermediate debriefs?
input('Have you scored the intermediate debriefs and added to wave_results?\n');

% Write in the wave_debriefing table
wave_debriefing.debrief_qc_pass = wave_results.debrief_qc_pass;

% Save all these files now:
save(fullfile(home,'results','analysis','waveData',...
    batch_wave_name,'wave_results_debrief_and_feedback.mat'),...
    'wave_debriefing','wave_int_feedback','wave_results');

fprintf('Wave feedback scoring saved!\n');

%% Now, update the columns in the long form data file

for iPtp = 1:height(wave_results)
    
    curr_ptp = wave_results.ptp{iPtp};
    
    % Find indices in the long_form_data
    idx_ptp = find(strcmp(wave_long_form_data.ptp,curr_ptp));
    
    % Update the feedback fields
    wave_long_form_data.debrief_qc_pass(idx_ptp) = wave_results.debrief_qc_pass(iPtp);
    wave_long_form_data.fb_int_qc_pass (idx_ptp) = wave_results.fb_int_qc_pass (iPtp);
    
end

% Save long form date
save(fullfile(home,'results','analysis','waveData',...
    batch_wave_name,'wave_long_form_data.mat'),'wave_long_form_data');

fprintf('Wave long_form_data updated with scores for debriefing and feedback.\n');