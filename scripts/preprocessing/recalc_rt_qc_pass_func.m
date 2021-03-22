function [long_form_data_all_ptp,results_table_all_ptp] = ...
    recalc_rt_qc_pass_func(saveCombinedData)

% Description:

% This function is called by preprocess.m but can also be used standalone. 
% It will take the data its passed or will load the datafiles for all ptps,
% and loop through each participant's data from each phase to see if their 
% mean RT is within 3SD of their group RT. It then saves these data files.

%% Global variables

if (nargin < 1)
    saveCombinedData = 1;
end

% Get the current directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end

% Load the files
[long_form_data_all_ptp,results_table_all_ptp,...
    debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
    load_all_participant_files(home);

%% Get only those participants that passed QCs up till now
long_form_data_qc_pass_ptp = get_only_qc_pass(long_form_data_all_ptp);

all_ptps = unique(long_form_data_qc_pass_ptp.ptp);

%% Create a table with the RT data for each congruency group
rt_table = table;
row_counter = 1;

for iCong = 0:1
    
    for iPhase = 1:2
        
        % Find all participants with this data
        iData_idx = long_form_data_qc_pass_ptp.congruency == iCong & ...
            long_form_data_qc_pass_ptp.phase == iPhase & ...
            long_form_data_qc_pass_ptp.session ~= -1;
        iData = long_form_data_qc_pass_ptp(iData_idx,:);
        
        % Get the RT distribution stats        
        compare_rt_mean  = nanmean(iData.rt);
        compare_rt_stdev = nanstd(iData.rt);
        
        lEdge = compare_rt_mean - 3 * compare_rt_stdev;
        uEdge = compare_rt_mean + 3 * compare_rt_stdev;      

        % Record in the rt table
        rt_table.congruency(row_counter) = iCong;
        rt_table.phase(row_counter)      = iPhase;
        
        rt_table.rt_mean(row_counter)    = compare_rt_mean;
        rt_table.rt_stdev(row_counter)    = compare_rt_stdev;
        rt_table.rt_lEdge(row_counter)    = lEdge;
        rt_table.rt_uEdge(row_counter)    = uEdge;        
        
        row_counter = row_counter + 1;
    end
    
    
end


%% Start the loop
for iPtp = 1:numel(all_ptps)
    curr_ptp = all_ptps{iPtp};
    
    for iiPhase = 1:2
        
        % Field variable to upate in the qc table and long_form_data
        iUpField = ['phase_' int2str(iiPhase) '_rt_qc_pass'];
        
        curr_ptp_phase_data_idx = strcmp(long_form_data_all_ptp.ptp,curr_ptp) & ...
            long_form_data_all_ptp.phase == iiPhase & ...
            long_form_data_all_ptp.session ~= -1;
        curr_ptp_phase_data     = long_form_data_all_ptp(curr_ptp_phase_data_idx,:);
        
        curr_ptp_phase_mean_rt  = nanmean(curr_ptp_phase_data.rt);
        
        curr_congruency   = curr_ptp_phase_data.congruency(1);
        
        % Get the RT distribution stats
        row_idx = find(rt_table.congruency == curr_congruency & ...
            rt_table.phase == iiPhase);

        lEdge = rt_table.rt_lEdge(row_idx);
        uEdge = rt_table.rt_uEdge(row_idx);
        
        %% Pass or fail?
        if curr_ptp_phase_mean_rt > uEdge || ...
                curr_ptp_phase_mean_rt < lEdge
            
            % Update the qc table
            results_table_all_ptp.(iUpField)(strcmp(results_table_all_ptp.ptp,curr_ptp))...
                = 0;
            
            % Update the long_form_data_all_ptp
            long_form_data_all_ptp.(iUpField)(strcmp(long_form_data_all_ptp.ptp,curr_ptp)) ...
                = zeros(nnz(strcmp(long_form_data_all_ptp.ptp,curr_ptp)),1);
            
            % Print a message
            fprintf(['Ptp ' curr_ptp ' RT of ' ...
                num2str(round(curr_ptp_phase_mean_rt)) ...
                ' is 3SD away from the group mean of ' ...
                num2str(round(compare_rt_mean)) ...
                '. \n3SD range is ' ...
                num2str(round(lEdge)) ' - ' ...
                num2str(round(uEdge))...
                '\n']);
        else
            % Update the results table
            results_table_all_ptp.(iUpField)(strcmp(results_table_all_ptp.ptp,curr_ptp))...
                = 1;
            
            % Update the long_form_data_all_ptp
            long_form_data_all_ptp.(iUpField)(strcmp(long_form_data_all_ptp.ptp,curr_ptp)) ...
                = ones(nnz(strcmp(long_form_data_all_ptp.ptp,curr_ptp)),1);
        end
    end % iPhase
end % iPtp

if saveCombinedData
    save(fullfile(home,'results','analysis',...
        'long_form_data_all_ptp.mat'),'long_form_data_all_ptp');
    save(fullfile(home,'results','analysis',...
        'results_table_all_ptp.mat'),'results_table_all_ptp');
end

end