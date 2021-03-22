%% Description

% Script to combine different JATOS files into one. Sometimes, I'd have to
% download them separately, because the JATOS website would crash
% otherwise.

%% General setup
dbstop if error;
close all;
clear; clc;

% Check you're in the right directory
home = pwd;

[~,name,~] = fileparts(home);

if ~strcmp('concept_gen_analysis',name)
    error('please change working directory to ./con_learn/unified_space/');
end

addpath(genpath(home));

saveLoc = fullfile(home,'data');

saveFileName = 'jatos_results_batch_3_wave_10.txt';

% sourceLoc = 'C:\Users\levan\Desktop\';
sourceLoc = fullfile(home,...
    'data','jatos_individual_exports','batch_3_wave_10');

%% Flags

saveFile = 0;


%% Specify file names

name_1 = fullfile(sourceLoc,'Study00426_ResultID-32283_WorkerID-15652_AUTO_EXPORT.jatos');
name_2 = fullfile(sourceLoc,'Study00426_ResultID-32284_WorkerID-15653_AUTO_EXPORT.jatos');
% name_3 = fullfile(sourceLoc,'jatos_results_20210127185411.txt');

%% Combine

% Read them
read_1 = readtable(name_1,'FileType','text');
read_2 = readtable(name_2,'FileType','text');
% read_2 = fileread(name_2);
% read_3 = fileread(name_3);
% 
% % If they are individual auto exports, shave off the initial text
% if strcmp(read_1(1:16),'Component Result')
%    
%     [file_split,~] = split(read_1,'[get_pid_comp_start---');
% 
%     
% end
% 
% 
% 
% combined_file = strcat(read_1,read_2);


%% Save

if saveFile
    
    save(fullfile(saveLoc,saveFileName),'combined_file');
    
end