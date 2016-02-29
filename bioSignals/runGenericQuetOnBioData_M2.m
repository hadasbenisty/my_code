% close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../gen_utils'));
addpath(genpath('../tiling'));

%% Init params
datapth = '..\..\..\datasets\biomed\M2';
files = { '4_4_14' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

rng(73631);

dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
%% Load Data


nt = 120;
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);


[nr, nt, nT] = size(X);
data = X;

%% Run Qu. 3D
params  = SetGenericQuestParamsM2;
params.row_emd.beta=.5;
params.col_emd.beta=.5;
params.trial_emd.beta=.5;
params.col_tree.treeDepth = 7;
params.row_tree.treeDepth = 7;
params.trial_tree.treeDepth = 7;
params.trial_tree.k = 2;
[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, permute(data, [2 3 1]) );
figure;
subplot(1,2,1);
[vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(1,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 100) 100*ones(1, 20) 200*ones(1, 240)], 'Time Colored by Tone');

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


[X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, nt);
[X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, nt);
[X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, nt);

tree_level = 2;
[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(row_tree{tree_level}, X, NeuronsLabels, NeuronsLabels);
[meanMat1, allMat1, meanMatAlltrials1] = getCentroidsByTree(row_tree{tree_level}, X1, NeuronsLabels, NeuronsLabels1);
[meanMat2, allMat2, meanMatAlltrials2] = getCentroidsByTree(row_tree{tree_level}, X2, NeuronsLabels, NeuronsLabels2);
[meanMat3, allMat3, meanMatAlltrials3] = getCentroidsByTree(row_tree{tree_level}, X3, NeuronsLabels, NeuronsLabels3);
figure;
for T = 1:size(meanMatAlltrials, 3)
    subplot(2,2,1);    imagesc(meanMatAlltrials(:,:,T));
    subplot(2,2,2);    imagesc(meanMatAlltrials1(:,:,T));
    subplot(2,2,3);    imagesc(meanMatAlltrials2(:,:,T));
    subplot(2,2,4);    imagesc(meanMatAlltrials3(:,:,T));
    pause;
end
[ aff_mat ] = CalcEmdAffOnTreeLevels( meanMat.', col_tree, tree_level, params.row_emd);
[vecs, vals] = CalcEigs(threshold(aff_mat, 0.0), 3);
[ row_order ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
% showing that the initial metric is  as good as the EMD
% 
[ aff_mat1,  ] = CalcInitAff( meanMat.', params.init_aff_row );
[vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);% 
[ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
plotByClustering(meanMat(row_order1, :), allMat(row_order1), ['8/12/14 Organized By Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat1(row_order1, :), allMat1(row_order1), ['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
plotByClustering(meanMat2(row_order1, :), allMat2(row_order1), ['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
plotByClustering(meanMat3(row_order1, :), allMat3(row_order1), ['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);


% 
% 
%  row_init_aff = feval(params.init_aff_row.initAffineFun, permute(data, [2 1 3]), params.init_aff_row);
% [row_tree_initial, row_embedding] = feval(params.row_tree.buildTreeFun, row_init_aff, params.row_tree);
% figure;plotTreeWithColors(row_tree_initial, 1:length(row_embedding));    title('Initial Nuerons Tree')
% tree_level = 2;
% [meanMat, allMat] = getCentroidsByTree(row_tree_initial{tree_level}, X, NeuronsLabels, NeuronsLabels);
% [meanMat1, allMat1] = getCentroidsByTree(row_tree_initial{tree_level}, X1, NeuronsLabels, NeuronsLabels1);
% [meanMat2, allMat2] = getCentroidsByTree(row_tree_initial{tree_level}, X2, NeuronsLabels, NeuronsLabels2);
% [meanMat3, allMat3] = getCentroidsByTree(row_tree_initial{tree_level}, X3, NeuronsLabels, NeuronsLabels3);
% [ aff_mat1,  ] = CalcInitAff( meanMat.', params.init_aff_row );
% [vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);% 
% [ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% % figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
% plotByClustering(meanMat(row_order1, :), allMat(row_order1), ['8/12/14 Organized By Tree Level ' num2str(tree_level)]);
% plotByClustering(meanMat1(row_order1, :), allMat1(row_order1), ['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
% plotByClustering(meanMat2(row_order1, :), allMat2(row_order1), ['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
% plotByClustering(meanMat3(row_order1, :), allMat3(row_order1), ['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
% 
