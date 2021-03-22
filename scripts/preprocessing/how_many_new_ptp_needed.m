function [long_form_data_all_ptp,results_table_all_ptp,...
            G,IID,...
            results_table_qc_pass,table_for_groups] = ...
            how_many_new_ptp_needed()

% Description:

%% Global variables

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
results_table_qc_pass = get_only_qc_pass(results_table_all_ptp);

all_ptps = unique(results_table_qc_pass.ptp);

%% Get all the combination of groups

% table to analyze
table_for_groups = results_table_qc_pass(:,{'congruency','concept_phase_1',...
    'concept_phase_2','arr_phase_1','arr_phase_2'});

% Turn the arrangements to a string, so findgroups can classify them
table_for_groups.arr_phase_1 = cellfun(@transpose,...
    table_for_groups.arr_phase_1,'UniformOutput',0);
table_for_groups.arr_phase_1 = cellfun(@num2str,...
    table_for_groups.arr_phase_1,'UniformOutput',0);
table_for_groups.arr_phase_2 = cellfun(@transpose,...
    table_for_groups.arr_phase_2,'UniformOutput',0);
table_for_groups.arr_phase_2 = cellfun(@num2str,...
    table_for_groups.arr_phase_2,'UniformOutput',0);

[G,IID] = findgroups(table_for_groups);

% Add the counts to the table
for iRow = 1:height(IID)
    IID.count(iRow) = sum(G == iRow);
end

end