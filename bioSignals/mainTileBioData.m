close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));

overwrite = false;
files = {'8_15_13_1-35' '8_12_14_1-40'};% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'
figspath1 = fullfile('D:\workWithBoss\summaries\D30',files{1});
for n=2:length(files)
    figspath1 = [figspath1 files{n}];
end
wrkspname = fullfile(figspath1, 'wrkspace.mat');
if exist(wrkspname, 'file') && ~overwrite
    load(wrkspname);
else
    rng(73631);
    %% Init params
    dorandperm_col = false;
    dorandperm_row = false;
    dorandperm_trials = false;
    eigsnum_col = 12;
    eigsnum_row = 12;
    eigsnum_trial= 12;
    row_alpha = .2;
    row_beta = 0;
    col_alpha = .2;
    col_beta = 0;
    trial_alpha = .2;
    trial_beta = 0;
    params  = SetQuest3DParams(eigsnum_col, eigsnum_row, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
    % params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trial);
    %% Load Data
    
    datapth = '..\..\..\datasets\biomed\D30';
    
    nt = 360;
    [X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);
    
    
    [nr, nt, nT] = size(X);
    % perm data
    if dorandperm_col
        col_perm = randperm(nt);
    else
        col_perm = 1:nt;
    end
    if dorandperm_row
        row_perm = randperm(nr);
    else
        row_perm = 1:nr;
    end
    if dorandperm_trials
        trial_perm = randperm(nT);
    else
        trial_perm = 1:nT;
    end
    data = X(row_perm, :, :);
    data = data(:, col_perm, :);
    data = data(:, :, trial_perm);
    
    
    %% Run Qu. 3D
    [row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data);
   
    %% Visualization
    mkNewFolder(figspath1);
    save(wrkspname);
end
% plotSlices(data, 9, 'Neurons','Time', 'Trials')
row_thresh = 0.0;
col_thresh = 0.0;
trials_thresh = 0;% 0.4
tonetimeS = 100;
tonetimeE = 120;
toneLabel = [ones(tonetimeS, 1); 2*ones(tonetimeE-tonetimeS, 1); 3*ones(size(col_aff, 1)-(tonetimeE), 1)];
savefigs=true;
if length(unique(expLabel))==1
    expLabel=[];
end

eigsnum_row = 3;
eigsnum_col = 3;
eigsnum_trials = 3;
% trim the trees so we'll have a much simple problem - just for debugging
trim_row_tree = row_tree(end-3:end);
trim_col_tree = col_tree(end-2:end);
trim_trials_tree = trials_tree(end-2:end);

volume = 7;
ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(data, 1), size(data, 2), 1);
tiling.isLeader = zeros(size(data, 1), size(data, 2), 1);
recursiveTiling2D(data(:,:,1), row_tree, col_tree, volume, ind2data, tiling, minErr);


% evaluate all trims having constant volume


figure;subplot(3,1,1);
treeplot(nodes(trim_row_tree),'.');title('Row Tree');
subplot(3,1,2);
treeplot(nodes(trim_col_tree),'.');title('Col Tree');
subplot(3,1,3);
treeplot(nodes(trim_trials_tree),'.');title('Trials Tree');