close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));

overwrite = false;
files = {'8_12_14_1-40'};% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

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
    
    mkNewFolder(figspath1);
    savefigs = true;
    save(wrkspname);
end
row_thresh = 0.0;
col_thresh = 0.0;
trials_thresh = 0;% 0.4
eigsnum_row = 3;
eigsnum_col = 3;
eigsnum_trials = 3;
tonetimeS = 100;
tonetimeE = 120;
toneLabel = [ones(tonetimeS, 1); 2*ones(tonetimeE-tonetimeS, 1); 3*ones(size(col_aff, 1)-(tonetimeE), 1)];
savefigs=true;
if length(unique(expLabel))==1
    expLabel=[];
end
% plotTreesAndEmbedding3D(figspath1, savefigs, 'Neurons', 'Time', 'Trials', ...
%     eigsnum_row, eigsnum_col, eigsnum_trials, ...
%     row_aff, col_aff, trials_aff, ...
%     row_thresh, col_thresh, trials_thresh, ...
%     row_tree, col_tree, trials_tree, ...
%     row_perm, col_perm, trial_perm,  [], toneLabel, expLabel);
figure;
subplot(1,2,1);
[vecs, vals] = CalcEigs(threshold(col_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(1,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 100) 100*ones(1, 20) 200*ones(1, 240)], 'Time Colored by Tone');

% figure;
% [vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
% figure;
% [vecs, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 4);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');
% subplot(3,2,6);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 35) 100*ones(1, 40) 150*ones(1, 45) 200*ones(1, 35) ], 'Trials Colored By Experiment');

figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trials_tree, 1:length(trials_aff));    title('Trialsl Tree');


[X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, nt);
[X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, nt);
[X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, nt);
tree_level = 3;
[meanMat, allMat] = getCentroidsByTree(row_tree{tree_level}, X, NeuronsLabels, NeuronsLabels);
[meanMat1, allMat1] = getCentroidsByTree(row_tree{tree_level}, X1, NeuronsLabels, NeuronsLabels1);
[meanMat2, allMat2] = getCentroidsByTree(row_tree{tree_level}, X2, NeuronsLabels, NeuronsLabels2);
[meanMat3, allMat3] = getCentroidsByTree(row_tree{tree_level}, X3, NeuronsLabels, NeuronsLabels3);


[ aff_mat ] = CalcEmdAffOnTreeLevels( meanMat.', col_tree, tree_level, params.row_emd);
[vecs, vals] = CalcEigs(threshold(aff_mat, 0.0), 3);
[ row_order ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:16, 'Nuerons Embedding');
% showing that the initial metric is  as good as the EMD
% 
[ aff_mat1,  ] = CalcInitAff( meanMat.', params.init_aff_row );
[vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);% 
[ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:16, 'Nuerons Embedding');
plotByClustering(meanMat(row_order1, :), allMat(row_order1), ['8/12/14 Organized By Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat1(row_order1, :), allMat1(row_order1), ['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
plotByClustering(meanMat2(row_order1, :), allMat2(row_order1), ['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
plotByClustering(meanMat3(row_order1, :), allMat3(row_order1), ['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
NeuronsLabels1(find(row_tree{tree_level}.folder_sizes==1))
