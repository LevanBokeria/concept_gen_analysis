%% Description

% Date: 20 April, 2021

% This script will take the real prolific IDs from the metadata file
% downloaded from prolific.co, and for each it will try to find a
% corresponding participant in the list of entered prolific IDs.

% The reason for making this script was that some participants messed up
% how they copy-pasted their prolific IDs, so their IDs in my raw data dont
% match their real prolific IDs. This created some issues in the pipeline,
% when trying to integrated prolific metadata with preprocessed and
% analysed data.

% This script is written after all the data was preprocessed, so its a
% post-hoc script. A new script should be written if we are going to get
% more data, such that preprocessing pipeline includes this mapping step.

%% General setup
clear; close all;
dbstop if error;

% Check you're in the right directory 
home = pwd;
[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/concept_gen_analysis/');
end
addpath(genpath(home));

%% Load the prolific meta data and the entered Prolific IDs from our preprocessed file

prolific_meta_data = readtable(fullfile(home,'data','prolific_meta_data',...
    'united_meta_data.xlsx'));

load(fullfile(home,'results','analysis','results_table_all_ptp.mat'));