%% Draft script to analyze study 2 data

% This script will try to fit a learning curve to each participants data.
% It will load the long form data table of all participant data and 
% use fminsearch to find the best fit of the model to data, using
% least squares error.

% It will record the output in the results_table_all_ptp.mat

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

%% Data parameters 
saveFiles          = 1;
drawNormalityPlots = 0;
plotFMSEstimation  = 0;

%% Load the data for all ptps:
[long_form_data_all_ptp,results_table_all_ptp,...
    debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
    load_all_participant_files(home);

nPtp = height(results_table_all_ptp);

%% Now, for each ptp, for each phase, fit the model.

%% Start the loop

% Create the fields to hold the data
results_table_all_ptp.phase_1_learning_rate_exp = NaN(nPtp,1);
results_table_all_ptp.phase_2_learning_rate_exp = NaN(nPtp,1);
results_table_all_ptp.phase_1_learning_rate_int = NaN(nPtp,1);
results_table_all_ptp.phase_2_learning_rate_int = NaN(nPtp,1);
results_table_all_ptp.phase_1_learning_rate_sse = NaN(nPtp,1);
results_table_all_ptp.phase_2_learning_rate_sse = NaN(nPtp,1);

for iPtp = 1:nPtp
    
    % Name of current participant
%     iPtp
    curr_ptp = results_table_all_ptp.ptp{iPtp};
    
    curr_congruency = results_table_all_ptp.congruency(iPtp);
    
    for iPhase = 1:2
        
        %% Get data for this ptp for this phase, minus practice trials
        idx_ptp      = strcmp(long_form_data_all_ptp.ptp,curr_ptp);
        idx_phase    = long_form_data_all_ptp.phase == iPhase;
        % Find practice trials
        idx_practice = long_form_data_all_ptp.session == -1;
        idx_wanted   = idx_ptp & idx_phase & ~idx_practice;
        
        % Some participants failed QC, so had no phase 2. Skip the loop for
        % them.
        if nnz(idx_wanted) == 0
            continue
        end
        
        curr_data    = long_form_data_all_ptp(idx_wanted,:);
        
        % Get the actual performance, correct/incorrect on each trial.
        % NaN means missed trials.
        y            = curr_data.correct;
        n_trials     = length(y);
            
        %% Calculate the intercept and learning rate using fminsearch
        
        % Starting value for c0
        c0         = 0.1;
        intercept0 = 0.5;
        % Optimize with respect to c
        params = c0;
        [fms_est_c_fixed_intercept,...
            fms_est_c_fixed_intercept_sse] = est_learning_rate(y,params,plotFMSEstimation);
        
        % Record everything in a table
        
        results_table_all_ptp.(...
            ['phase_' int2str(iPhase) '_learning_rate_exp'])(iPtp) = ...
            fms_est_c_fixed_intercept;
        results_table_all_ptp.(...
            ['phase_' int2str(iPhase) '_learning_rate_int'])(iPtp) = ...
            intercept0;
        results_table_all_ptp.(...
            ['phase_' int2str(iPhase) '_learning_rate_sse'])(iPtp) = ...
            fms_est_c_fixed_intercept_sse;        
    end % Phase
end % iPtp

%% Now, calculate differences between phase 1 and phase 2 for each dep var

% For learning rate
results_table_all_ptp.phase_2_min_phase_1_learning_rate_exp = ...
    results_table_all_ptp.phase_2_learning_rate_exp - ...
    results_table_all_ptp.phase_1_learning_rate_exp;

% For ses 1 and 2 performance
results_table_all_ptp.phase_2_min_phase_1_ses_1_2_perf = ...
    results_table_all_ptp.phase_2_ses_1_2_perf - ...
    results_table_all_ptp.phase_1_ses_1_2_perf;

%% Get only qc pass ptps
results_table_qc_pass_ptp = get_only_qc_pass(results_table_all_ptp);

%% Define the dependent variables as separate variables
% Learning rate
learning_rate_congruent_0 = results_table_qc_pass_ptp.phase_2_min_phase_1_learning_rate_exp(...
    results_table_qc_pass_ptp.congruency == 0);
learning_rate_congruent_1 = results_table_qc_pass_ptp.phase_2_min_phase_1_learning_rate_exp(...
    results_table_qc_pass_ptp.congruency == 1);

% Variances 


% ses 1-2 performance
ses_1_2_perf_congruent_0 = results_table_qc_pass_ptp.phase_2_min_phase_1_ses_1_2_perf(...
    results_table_qc_pass_ptp.congruency == 0);
ses_1_2_perf_congruent_1 = results_table_qc_pass_ptp.phase_2_min_phase_1_ses_1_2_perf(...
    results_table_qc_pass_ptp.congruency == 1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transform the data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add a constant
learning_rate_congruent_0_positive = learning_rate_congruent_0 + 22;
learning_rate_congruent_1_positive = learning_rate_congruent_1 + 22;

%% Take the log
learning_rate_congruent_0_log = log(learning_rate_congruent_0_positive);
learning_rate_congruent_1_log = log(learning_rate_congruent_1_positive);

%% Square root
learning_rate_congruent_0_sqrt = sqrt(learning_rate_congruent_0_positive);
learning_rate_congruent_1_sqrt = sqrt(learning_rate_congruent_1_positive);

%% Power law 3
learning_rate_congruent_0_cube_root = ...
    learning_rate_congruent_0_positive.^(1/3);
learning_rate_congruent_1_cube_root = ...
    learning_rate_congruent_1_positive.^(1/3);

%% Power law 4
learning_rate_congruent_0_four_root = ...
    learning_rate_congruent_0_positive.^(1/4);
learning_rate_congruent_1_four_root = ...
    learning_rate_congruent_1_positive.^(1/4);

%% Reciprocal
learning_rate_congruent_0_recip = learning_rate_congruent_0.^(-1);
learning_rate_congruent_1_recip = learning_rate_congruent_1.^(-1);


%% %%%%%%%%%%%%%%%%%%%%%% Plots for normality %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if drawNormalityPlots

    %% Learning rate:
    histBins = 25;
    xMin = -30;
    xMax = 30;
    yMin = 0;
    yMax = 50;
    titleText1 = 'Incongruent';
    titleText2 = 'Congruent';
    supTitleText = 'Learning Rate';

    % Plot histograms
    drawHist(learning_rate_congruent_0,...
        learning_rate_congruent_1,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Log Learning Rate';
    drawHist(learning_rate_congruent_0_log,...
        learning_rate_congruent_1_log,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Sqrt Learning Rate';
    drawHist(learning_rate_congruent_0_sqrt,...
        learning_rate_congruent_1_sqrt,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Cube Root Learning Rate';
    drawHist(learning_rate_congruent_0_cube_root,...
        learning_rate_congruent_1_cube_root,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Four Root Learning Rate';
    drawHist(learning_rate_congruent_0_four_root,...
        learning_rate_congruent_1_four_root,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Reciprocal Learning Rate';
    drawHist(learning_rate_congruent_0_recip,...
        learning_rate_congruent_1_recip,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    % Plot QQ plots
    supTitleText = 'Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0,...
        learning_rate_congruent_1,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Log Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0_log,...
        learning_rate_congruent_1_log,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Sqrt Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0_sqrt,...
        learning_rate_congruent_1_sqrt,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Cube Root Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0_cube_root,...
        learning_rate_congruent_1_cube_root,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Four Root Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0_four_root,...
        learning_rate_congruent_1_four_root,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);

    supTitleText = 'Reciprocal Learning Rate QQ Plots';
    drawQQPlot(learning_rate_congruent_0_recip,...
        learning_rate_congruent_1_recip,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText);
    %% Session 1-2 performance:
    supTitleText = 'Ses 1-2 avg performance';

    drawHist(ses_1_2_perf_congruent_0,...
        ses_1_2_perf_congruent_1,...
        histBins,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText)

    % Plot QQ plots
    supTitleText = 'QQ Plots';

    drawQQPlot(ses_1_2_perf_congruent_0,...
        ses_1_2_perf_congruent_1,...
        xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText)

end
%% Save the data
if saveFiles
    save(fullfile(home,'results','analysis','results_table_qc_pass_ptp_analyzed.mat'),...
        'results_table_qc_pass_ptp');
    
    % Save as excel
    results_table_qc_pass_ptp_excel = results_table_qc_pass_ptp;
    
    results_table_qc_pass_ptp_excel.targetNamesUsed_phase_1 = [];
    results_table_qc_pass_ptp_excel.targetNamesUsed_phase_2 = [];
    
    writetable(results_table_qc_pass_ptp_excel,...
        fullfile(home,'results','analysis','results_table_qc_pass_ptp_analyzed.xlsx'));
    
end


function drawHist(data1, data2, ...
    histBins,xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText)

figure
subplot(1,2,1)

histogram(data1);

% xlim([xMin xMax]);
% ylim([yMin yMax]);
title(titleText1)

subplot(1,2,2)

histogram(data2);

% xlim([xMin xMax]);
% ylim([yMin yMax]);
title(titleText2)

suptitle(supTitleText)

end

function drawQQPlot(data1, data2, ...
    xMin,xMax,yMin,yMax,titleText1,titleText2,supTitleText)

figure
subplot(1,2,1)

qqplot(data1);

% xlim([xMin xMax]);
% ylim([yMin yMax]);
title(titleText1)

subplot(1,2,2)

qqplot(data2);

% xlim([xMin xMax]);
% ylim([yMin yMax]);
title(titleText2)

suptitle(supTitleText)

end