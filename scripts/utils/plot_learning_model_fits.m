%% Description

% This script just plot the data and model fit of specific participants
% that had a really large exponent estimate. The purpose is to see why that
% is the case.

%% Global parameters
clear; close all;

dbstop if error;

% Check you're in the right directory 
home = pwd;
   
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

addpath(genpath(home));

%% Load the data for all ptps:
[long_form_data_all_ptp,results_table_all_ptp,...
    debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
    load_all_participant_files(home);

results_table_qc_pass_ptp_analyzed = load(fullfile(home,'results',...
    'analysis','results_table_qc_pass_ptp_analyzed.mat'));
results_table_qc_pass_ptp_analyzed = results_table_qc_pass_ptp_analyzed.results_table_qc_pass_ptp;

%% Find the ones with large estimates
idx_large_est = results_table_qc_pass_ptp_analyzed.phase_1_learning_rate_exp > 1 | ...
    results_table_qc_pass_ptp_analyzed.phase_2_learning_rate_exp > 1;

results_table_large_est_ptp = results_table_qc_pass_ptp_analyzed(idx_large_est,:);

nPtp = height(results_table_large_est_ptp);
%% Start the loop for each participant
for iPtp = 1:nPtp
    
    curr_ptp = results_table_large_est_ptp.ptp{iPtp};
    
    curr_est_p1 = results_table_large_est_ptp.phase_1_learning_rate_exp(iPtp);
    curr_est_p2 = results_table_large_est_ptp.phase_2_learning_rate_exp(iPtp);

    % Find the row in long form data
    idx_ptp_data_p1 = find(strcmp(long_form_data_all_ptp.ptp,curr_ptp) & ...
        long_form_data_all_ptp.phase == 1);
    idx_ptp_data_p2 = find(strcmp(long_form_data_all_ptp.ptp,curr_ptp) & ...
        long_form_data_all_ptp.phase == 2);
    
    curr_data_p1 = long_form_data_all_ptp.correct(idx_ptp_data_p1);
    curr_data_p2 = long_form_data_all_ptp.correct(idx_ptp_data_p2);
% 
    trials = 1:length(curr_data_p1);
    model_fit_data_p1 = 1 - (0.5 * exp(-curr_est_p1 .* (trials-1)));
    trials = 1:length(curr_data_p2);
    model_fit_data_p2 = 1 - (0.5 * exp(-curr_est_p2 .* (trials-1)));

    % Draw
    drawModelFit(curr_data_p1,curr_data_p2,...
        model_fit_data_p1,model_fit_data_p2,curr_ptp);
    
    
end


%% Now, for each ptp, for each phase, fit the model.
% Define the figures
fig_fit_plots = figure(2);
% set(fig_fit_plots,'Position',[548   117   794   840]);
set(fig_fit_plots,'Position',[274    72   794   840]);


function drawModelFit(curr_data_p1,curr_data_p2,...
    model_fit_data_p1,model_fit_data_p2,curr_ptp)

figure
set(gcf,'Position',[-175  -488   794   241]);

subplot(1,2,1);
plot(1:length(curr_data_p1),curr_data_p1,'LineWidth',1);
hold on;
plot(1:length(curr_data_p1),model_fit_data_p1,'LineWidth',1.5);

grid on;
xlim([0 170]);
ylim([-0.1 1.1]);
xlabel('trials');
ylabel('Accuracy');
title('Phase 1');

subplot(1,2,2);
plot(1:length(curr_data_p2),curr_data_p2,'LineWidth',1);
hold on;
plot(1:length(curr_data_p2),model_fit_data_p2,'LineWidth',1.5);

grid on;
xlim([0 170]);
ylim([-0.1 1.1]);
xlabel('trials');
ylabel('Accuracy');
title('Phase 2');


suptitle(curr_ptp);


end

