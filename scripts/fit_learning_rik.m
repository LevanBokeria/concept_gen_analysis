%% Draft script to analyze study 2 data

% This script will try to fit a learning curve to each participants data.
% It will first loop through all participants given to it, to preprocess their
% data. Creates a table of neatly organized performance data.
% Then, it uses fminsearch to find the best fit of the model to data, using
% least squares error

% Flags allow choosing to plot all participants together, or in separate
% figures. Also, whether to plot all toys data or specific toys.
% Additional parameters to choose how to smooth data

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

% Variables 
% - Participants: data for each participant separately in a different cell.
% - data_long: long-format table containing all of data from all participants in one big table. 
load('participant_data.mat');

nPtp = height(participants);

savePerfFigure = 0;
averageAcrossPtp = 0;

% Flag to draw plots
drawPlots         = 0;
plotFMSEstimation = 0;
writeExcel        = 0;

%% Now, for each ptp, for each phase, fit the model.

if drawPlots
    % Define the figures
    fig_c_search  = figure(1);
    fig_fit_plots = figure(2);
    set(fig_c_search,'Position',[548   117   794   840]);
    set(fig_fit_plots,'Position',[548   117   794   840]);
end
%% Create predicted dataset for various parameters 

% intercept for the learning model
intercept_fixed = 0.5;

% Define range of exponent values
c_vector = 0:0.001:20;

% Create a matrix of nTrials - by - c, containing expected values for each
% trial varied by the exponent value
% nTrials set to 200 for no reason, just has to be at least 4*42.
y_hat = [];
for iT = 1:300
    y_hat(iT,:) = 1 - intercept_fixed * exp(-c_vector*(iT-1));
end

% Create a table to hold the analyzed data.
fit_table = table;

% Just a counter to write in the table
ctr = 1; 

%% Start the loop

for iPtp = 21:nPtp
    
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
        
        % Get the actual performance, correct/incorrect on each trial.
        % NaN means missed trials.
        y            = curr_data.correct;
        n_trials     = length(y);
        
        %% Calculate least squares for each predicted performance based on c
        % value.
        
            % Create a metrix of the data replicated many times 
            y_rep = repmat(y,1,length(c_vector));
        
            % Truncate the y_hat to have as many trials as the real data
            y_hat_trunk = y_hat(1:n_trials,:);
            
            % Subtract from y and square the error
            y_hat_trunk_error_squared = (abs(y_rep - y_hat_trunk)).^2;
            
            % Sum of squared error. one for each c coefficient
            sse = nansum(y_hat_trunk_error_squared,1);
            
            % Find the best c value, with minimal sse
            [min_sse,idx_min_sse] = min(sse);
            best_c = c_vector(idx_min_sse);
            
        %% Calculate the intercept and learning rate using fminsearch
        
        % Starting value for c0
        c0         = 0.1;
        intercept0 = 0.5;
        
        % 1. try with fixed intercept 0.5
            
            % Optimize with respect to c
            params = c0;
            [fms_est_c_fixed_intercept,...
             fms_est_c_fixed_intercept_sse] = est_learning_rate(y,params,plotFMSEstimation);
            
        % 2. Also optimize both, just to compare
            params = [c0,intercept0];
            [out_params,fms_est_both_sse] = est_learning_rate(y,params,plotFMSEstimation);
            
            fms_est_both_c      = out_params(1);
            fms_est_both_intercept = out_params(2);
            
        %% Calculate number of trials needed to reach criterion on every toy
        
        % Reset variables
        rowidx_of_last_trial = [];
        n_trials_all_targets = [];
        
        % Using a trailing sliding window of size w and criterion crit
        w    = 10;
        crit = 0.85;
        
        % Add rowindex to the curr_data, to use in finding trial index
        % by which criterion was reached for every toy.
        curr_data.rowidx = (1:height(curr_data))';
        
        % Which prompt toys were used? Get their names
        [prompt_groups, prompt_names] = findgroups(curr_data.prompt_img_name);
        
        % For each toy, see when the criterion of 85% was reached.
        for iGroup = 1:numel(prompt_names)
            curr_prompt = prompt_names{iGroup};
            
            % Subselect just the data on that prompt toy
            iTable = curr_data(strcmp(curr_data.prompt_img_name,curr_prompt),:);
            
            % Now for this prompt, at which trial did they reach criterion
            n_trials_took = sliding_window(iTable.correct,w,crit);
            
            % Which trial was that in the whole experiment?
            if isnan(n_trials_took)
                % So for this toy, they never reached the criterion
                rowidx_of_last_trial(iGroup) = n_trials_took;
            else
                rowidx_of_last_trial(iGroup) = iTable.rowidx(n_trials_took);
            end
        end
        
        % By which trial did they reach crit for the last toy?
        if nnz(isnan(rowidx_of_last_trial)) > 0
            n_trials_all_targets = NaN;
        else
            n_trials_all_targets = max(rowidx_of_last_trial);
        end
        
    
            %% Plots
            if drawPlots
                % Plot the search for minimum c
                figure(fig_c_search)
                subplot(nPtp,2,ctr);
                plot(sse,'LineWidth',1.5);
                grid on
                xlim([-10,length(c_vector)]);
                ylim([0 50]);
                xlabel('c value');
                ylabel('LSE');
                title(curr_ptp);
                
                % Plot the data and the learning curve based on best c
                figure(fig_fit_plots)
                subplot(nPtp,2,ctr);
                plot(1:length(y),y,'LineWidth',1);
                hold on;
                plot(1:length(y),y_hat_trunk(:,idx_min_sse),'LineWidth',1.5);
                
                grid on;
                xlim([0 170]);
                ylim([-0.1 1.1]);
                xlabel('trials');
                ylabel('Accuracy');
                title(curr_ptp);
            end
                        
        %% Record everything in a table
        fit_table.ptp{ctr}               = curr_ptp;
        fit_table.congruency(ctr)        = participants.congruency(iPtp);
        fit_table.phase(ctr)             = iPhase;
        
        % Record brute force search results
        fit_table.bf_c(ctr)                 = best_c;
        fit_table.bf_sse(ctr)               = min_sse;
%         fit_table.bf_squared_error_mat{ctr} = y_hat_trunk_error_squared;
        fit_table.bf_sse_for_each_c{ctr}    = sse;
        fit_table.bf_c_vector{ctr}          = c_vector;
        
        % Record fminsearch results estimating both c and intercept
        fit_table.fms_est_both_c(ctr)         = fms_est_both_c;
        fit_table.fms_est_both_intercept(ctr) = fms_est_both_intercept;
        fit_table.fms_est_both_sse(ctr)       = fms_est_both_sse;
        
        % Record fminsearch with fixed intercept
        fit_table.fms_est_c_fixed_intercept(ctr)     = fms_est_c_fixed_intercept;
        fit_table.fms_fixed_intercept(ctr)           = intercept_fixed;
        fit_table.fms_est_c_fixed_intercept_sse(ctr) = fms_est_c_fixed_intercept_sse;
        
        % Record the c and intercept to use in final analysis
%         if iPhase == 1
%             fit_table.c_to_analyze(ctr)         = fms_est_c_fixed_intercept;
%             fit_table.intercept_to_analyze(ctr) = intercept_fixed;
%         elseif iPhase == 2
%             fit_table.c_to_analyze(ctr)         = fms_est_both_c;
%             fit_table.intercept_to_analyze(ctr) = fms_est_both_intercept;
%         end
        
        % ses 1 perf:
        curr_data_ses_1 = curr_data(curr_data.session == 1,:);
        fit_table.session_1_perf(ctr) = nanmean(curr_data_ses_1.correct);

        % ses 2 perf:
        curr_data_ses_2 = curr_data(curr_data.session == 2,:);
        fit_table.session_2_perf(ctr) = nanmean(curr_data_ses_2.correct);
        
        % ses 1 and 2 performance 
        fit_table.session_1_2_perf(ctr) = ...
            nanmean([curr_data_ses_1.correct; curr_data_ses_2.correct]);
        
        % Trials to criterion
        fit_table.n_trials_to_crit(ctr) = n_trials_all_targets;
        
        % Iterate the counter 
        ctr = ctr + 1;
    end
end

%% Delete empty rows in the fit_table
idx_empty = cellfun(@isempty,fit_table.ptp);
fit_table(idx_empty,:) = [];

%% Set figure titles
if drawPlots
    figure(fig_c_search)
    sgtitle('Finding the best exponent');
    figure(fig_fit_plots)
    sgtitle('Data and the fit');


%     figure,plot(fit_table.session_1_perf,fit_table.c_to_analyze,'o')
%     [R,p]=corr(fit_table.session_1_perf,fit_table.c_to_analyze)
    
    
%     figure,plot(fit_table.mean_p(fit_table.phase==1),fit_table.exp(fit_table.phase==1),'o')
%     [R,p]=corr(fit_table.mean_p(fit_table.phase==1),fit_table.exp(fit_table.phase==1))

end

if writeExcel
    
    fit_table_to_write = fit_table;
    fit_table_to_write.bf_sse_for_each_c = [];
    fit_table_to_write.bf_c_vector       = [];
    
    writetable(fit_table_to_write,'fit_table1.xlsx');
    
end
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