function [long_form_data_all_ptp,results_table_all_ptp,...
    debrief_tbl_all_ptp, fb_inter_tbl_all_ptp] = ...
    load_all_participant_files(home)

% Long form data
if exist(fullfile(home,'results','analysis',...
        'long_form_data_all_ptp.mat'))
    % Load it
    long_form_data_all_ptp = load(fullfile(home,'results',...
        'analysis','long_form_data_all_ptp.mat'));
    long_form_data_all_ptp = long_form_data_all_ptp.long_form_data_all_ptp;
    
    fprintf('long_form_data_all_ptp already existed. Loading it...\n');
else
    % Make a fresh one
    long_form_data_all_ptp = table;
    fprintf('long_form_data_all_ptp created afresh!\n');
end

% QC table
if exist(fullfile(home,'results','analysis',...
        'results_table_all_ptp.mat'))
    % Load it
    results_table_all_ptp = load(fullfile(home,'results',...
        'analysis','results_table_all_ptp.mat'));
    results_table_all_ptp = results_table_all_ptp.results_table_all_ptp;
    
    fprintf('results_table_all_ptp already existed. Loading it...\n');
else
    % Make a fresh one
    results_table_all_ptp = table;
    fprintf('results_table_all_ptp created afresh!\n');
end

% Final debriefs file
if exist(fullfile(home,'results','analysis',...
        'debrief_tbl_all_ptp.mat'))
    % Load it
    debrief_tbl_all_ptp = load(fullfile(home,'results',...
        'analysis','debrief_tbl_all_ptp.mat'));
    debrief_tbl_all_ptp = debrief_tbl_all_ptp.debrief_tbl_all_ptp;
    
    fprintf('debrief_tbl_all_ptp already existed. Loading it...\n');
else
    % Make a fresh one, with field names predefined.
    debrief_tbl_all_ptp                              = table;
%     debrief_tbl_all_ptp.ptp                      {1} = {};
%     debrief_tbl_all_ptp.end_point                {1} = {};
%     debrief_tbl_all_ptp.debrief_qc_pass          (1) = NaN;
%     debrief_tbl_all_ptp.general_questions_page_Q0{1} = {};
%     debrief_tbl_all_ptp.general_questions_page_Q1{1} = {};
%     debrief_tbl_all_ptp.general_questions_page_Q2{1} = {};
%     debrief_tbl_all_ptp.general_questions_page_Q3{1} = {};
%     debrief_tbl_all_ptp.general_questions_page_Q4{1} = {};
%     debrief_tbl_all_ptp.phase_1_debriefer_Q0     {1} = {};
%     debrief_tbl_all_ptp.phase_1_debriefer_Q1     {1} = {};
%     debrief_tbl_all_ptp.phase_2_debriefer_Q0     {1} = {};
%     debrief_tbl_all_ptp.phase_2_debriefer_Q1     {1} = {};
%     debrief_tbl_all_ptp.phase_1_spatial_layout_Q0{1} = {};
%     debrief_tbl_all_ptp.phase_2_spatial_layout_Q0{1} = {};
%     debrief_tbl_all_ptp.phase_3_was_1_page_Q0    {1} = {};
%     debrief_tbl_all_ptp.phase_3_was_1_page_Q1    {1} = {};
%     
    fprintf('debrief_tbl_all_ptp created afresh!\n');
end

% Intermediate debriefs file
if exist(fullfile(home,'results','analysis',...
        'fb_inter_tbl_all_ptp.mat'))
    % Load it
    fb_inter_tbl_all_ptp = load(fullfile(home,'results',...
        'analysis','fb_inter_tbl_all_ptp.mat'));
    fb_inter_tbl_all_ptp = fb_inter_tbl_all_ptp.fb_inter_tbl_all_ptp;
    
    fprintf('fb_inter_tbl_all_ptp already existed. Loading it...\n');
else
    % Make a fresh one
    fb_inter_tbl_all_ptp = table;
    
    fprintf('fb_inter_tbl_all_ptp created afresh!\n');
end
end