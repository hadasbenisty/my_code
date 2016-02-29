% close all;
% run generic trees using correlation clustering using Shai Bagon's code
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../gen_utils'));
addpath(genpath('../tiling'));
addpath(genpath('../../LargeScaleCC1.0/'));
files = {'8_6_14_1-20_control' '8_6_14_21-60_cno'};% '8_6_14_21-60_cno' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
%% Load Data

datapth = '..\..\..\datasets\biomed\D8';

nt = 120;
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);


[nr, nt, nT] = size(X);
data = X;

%% Run Qu. 3D
params  = SetGenericQuestParamsD30;
params.row_emd.beta=.5;
params.col_emd.beta=.5;
params.trial_emd.beta=.5;
params.col_tree.treeDepth = 7;
params.row_tree.treeDepth = 7;
params.trial_tree.treeDepth = 7;
params.trial_tree.k = 2;
params.col_tree.clusteringAlgo = @LargeScaleCorrWrapper;
params.row_tree.clusteringAlgo = @svdClassWrapper;
params.trial_tree.clusteringAlgo = @LargeScaleCorrWrapper;
params.col_tree.runOnEmbdding = true;
params.row_tree.runOnEmbdding = true;
params.trial_tree.runOnEmbdding = true;
params.init_aff_col.RangeMinus1to1 = false;
params.init_aff_row.RangeMinus1to1 = false;
params.init_aff_trial.RangeMinus1to1 = false;
params.col_tree.ccAlgo = 'expand_explore';
params.row_tree.ccAlgo = 'swap_explore';
params.trial_tree.ccAlgo = 'expand_explore';
[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, permute(data, [2 3 1]) );

figure;
subplot(3,2,1);
[vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(3,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');

subplot(3,1,2);
[vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
subplot(3,1,3);
[vecs, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');

figure;    subplot(3,1,1);
plotTreeWithColors(trial_tree, 1:length(trial_dual_aff));    title('Trials Tree');
subplot(3,1,2);    plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Time Tree')
subplot(3,1,3);    plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Nuerons Tree');


