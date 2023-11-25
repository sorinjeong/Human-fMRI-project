clear all; clc;
addpath(genpath('D:\leelab\Human fMRI projects\OCAT_DIR'));

%% 
file_path_in = '../data/data_bhv_raw';
file_path_out = '../data/data_bhv_result'; 

%% 
n_sbj = 31;

for sbj_i = 1: n_sbj
    c_sbj = num2str(sbj_i, '%02.f');







%% make output directory
   path_out = {};
   path_out{end+1} = fullfile(file_path_out,'individual',c_sbj);
   path_out{end+1} = fullfile(file_path_out,'total');
   path_out{end+1} = fullfile(file_path_out,'GLM');
   path_out{end+1} = fullfile(file_path_out,'figure');
     
    if ~exist(file_path_out,"dir")
      mkdir(file_path_out);
      for i=1:length(path_out)
         mkdir(path_out{i});
      end
    end

