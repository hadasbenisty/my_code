close all;
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
run2D = false;
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

%% Run Qu. 2D  - comparing the generic dims to 2D - both with BU

% params1BU = SetGenericDimsQuestParams(2, false);
% [ TreesBU, dual_affBU, init_affBU ] = RunGenericDimsQuestionnaire( params1BU, data(:,:,1) );
% figure;
% subplot(1,2,1);
% [vecs, vals] = CalcEigs(threshold(dual_affBU{2}, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
% subplot(1,2,2);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');
% 
% figure;
% [vecs, vals] = CalcEigs(threshold(dual_affBU{1}, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding')
% 
% 
% paramsBU  = SetQuestParams(0,0,0,0);
% paramsBU.init_aff_row.on_rows = false;
% paramsBU.init_aff_col = paramsBU.init_aff;
% paramsBU.col_tree.eigs_num = 50;
% paramsBU.row_tree.eigs_num = 50;
% [ row_treeBU, col_treeBU, row_dual_affBU, col_dual_affBU ] = RunQuestionnaire( paramsBU, data(:,:,1) );
% figure;
% subplot(1,2,1);
% [vecs, vals] = CalcEigs(threshold(col_dual_affBU, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
% subplot(1,2,2);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');
% 
% figure;
% [vecs, vals] = CalcEigs(threshold(row_dual_affBU, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding')
% 
% params1 = SetGenericDimsQuestParams(2, true);
% [ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params1, data(:,:,1) );
% figure;
% subplot(1,2,1);
% [vecs, vals] = CalcEigs(threshold(dual_aff{1}, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
% subplot(1,2,2);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');
% 
% figure;
% [vecs, vals] = CalcEigs(threshold(dual_aff{2}, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding')
% params  = SetGenericQuestParamsD30;
% params.init_aff_col.initAffineFun= @CalcInitAff;
% params.init_aff_row.initAffineFun = @CalcInitAff;
% params.init_aff_trial.initAffineFun = @CalcInitAff;
% params.col_tree.CalcAffFun =  @CalcEmdAff;
% params.row_tree.CalcAffFun =  @CalcEmdAff;
% [ col_tree, row_tree, col_dual_aff, row_dual_aff  ] = RunGenericQuestionnaire2D( params, permute(data(:,:,1), [2 1 3]) );
% figure;
% subplot(1,2,1);
% [vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
% subplot(1,2,2);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');
% 
% figure;
% [vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
% plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding')


%% Run Qu. 3D - comparing the generic dims to 3D
params1BU = SetGenericDimsQuestParams(3, false);
[ TreesBU, dual_affBU, init_affBU ] = RunGenericDimsQuestionnaire( params1BU, data );
 paramsBU  = SetQuest3DParams(50, 50, 50, 1, 0, 1, 0, 1, 0 );

[ row_treeBU, col_treeBU, trials_treeBU, row_dual_aff_finalBU, col_dual_aff_finalBU, trials_dual_aff_finalBU, col_init_affBU, row_init_affBU, trials_init_affBU ] = RunQuestionnaire3D( paramsBU, permute(data, [2 1 3]) );

figure;
subplot(1,2,1);
[vecs, vals] = CalcEigs(threshold(row_dual_aff_finalBU, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(1,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');

figure;
[vecs, vals] = CalcEigs(threshold(col_dual_aff_finalBU, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
figure;
[vecs, vals] = CalcEigs(threshold(trials_dual_aff_finalBU, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');

figure;    subplot(3,1,1);
plotTreeWithColors(col_treeBU, 1:length(col_dual_aff_finalBU));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_treeBU, 1:length(row_dual_aff_finalBU));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trials_treeBU, 1:length(trials_dual_aff_finalBU));    title('Trialsl Tree');



params1 = SetGenericDimsQuestParams(3, true);
[ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params1, data );
figure;    subplot(3,1,1);
plotTreeWithColors(Trees{1}, 1:length(col_dual_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(Trees{2}, 1:length(row_dual_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(Trees{3}, 1:length(trial_dual_aff));    title('Trials Tree');

params  = SetGenericQuestParamsD30;

params.col_tree.CalcAffFun =  @CalcEmdAff3D;
params.row_tree.CalcAffFun =  @CalcEmdAff3D;

params.init_aff_col.initAffineFun= @CalcInitAff3D;
params.init_aff_row.initAffineFun = @CalcInitAff3D;
params.init_aff_trial.initAffineFun = @CalcInitAff3D;
params.init_aff_trial.initAffineFun = @CalcInitAff3D;
[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, permute(data, [2 3 1]) );

figure;
subplot(1,2,1);
[vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(1,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');

figure;
[vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
figure;
[vecs, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');

figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trial_tree, 1:length(trial_dual_aff));    title('Trialsl Tree');

