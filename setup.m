%% Setup ReSpace-MOBO repository

clear;
clc;

% Get repository root
projectRoot = fileparts(mfilename('fullpath'));

% Add repository and subfolders to MATLAB path
addpath(genpath(projectRoot));

disp('ReSpace-MOBO repository setup complete.');