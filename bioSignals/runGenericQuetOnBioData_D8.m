close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));

files = {'8_6_14_1-20_control' '8_6_14_21-60_cno' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
eigsnum_col = 12;eigsnum_row = 12;eigsnum_trial= 12;
row_alpha = .2;row_beta = 0;col_alpha = .2;col_beta = 0;trial_alpha = .2;trial_beta = 0;
params  = SetQuest3DParams(eigsnum_row, eigsnum_col, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
%% Load Data

datapth = '..\..\..\datasets\biomed\D8';

nt = 120;
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);


[nr, nt, nT] = size(X);
data = X;

%% Run Qu. 2D
params  = SetGenericQuestParams;
   
[ row_tree, col_tree, trial_tree, row_dual_aff, col_dual_aff, trial_dual_aff ] = RunGenericQuestionnaire3D( params, data );
figure;
subplot(3,2,1);
[vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.2), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(3,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');

subplot(3,1,2);
[vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
subplot(3,1,3);
[vecs, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');
% subplot(3,2,6);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 20) 100*ones(1, 40)], 'Trials Colored by CNO');
figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trial_tree, 1:length(trial_dual_aff));    title('Trialsl Tree');

% folders = unique(row_tree{4}.clustering);  
% for ci = 1:length(folders)
%     inds2folder = find(row_tree{3}.clustering == folders(ci));
%     R = ceil(sqrt(length(inds2folder)));
%     figure;
%     for r = 1:length(inds2folder)
%        subplot(R, R, r);
%        imagesc(permute(data(inds2folder(r), :, :), [3 2 1]));
%        title(NeuronsLabels{inds2folder(r)});
%     end
% end
    
    
    % %% Order by tree
% 
% [~, row_order] = sort(row_tree{2}.clustering);
% [~, col_order] = sort(col_tree{2}.clustering);
% orderedData = data(row_order, :);
% orderedData = orderedData(:, col_order);
% 
% % prepare trees for recursion
% for treeLevel = 1:length(row_tree)
%     row_orderedtree{treeLevel} = row_tree{treeLevel};
%     row_orderedtree{treeLevel}.clustering = row_tree{treeLevel}.clustering( row_order);
% end
% for treeLevel = 1:length(col_tree)
%     col_orderedtree{treeLevel} = col_tree{treeLevel};
%     col_orderedtree{treeLevel}.clustering = col_tree{treeLevel}.clustering( col_order);
% end
% 
% figure;
% subplot(2,2,1);
% imagesc(data);title('Data');
% subplot(2,2,2);
% imagesc(orderedData);title('Ordered Data');
% subplot(2,2,3);
% plotTreeWithColors(row_orderedtree, 1:size(data,1));
% title('Row Tree');
% 
% subplot(2,2,4);
% plotTreeWithColors(col_orderedtree, 1:size(data,2))
% title('Col Tree');
% 
% 
% ind2data = 1;
% minErr = Inf;
% tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
% tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
% solutionTiling = [];
% figure;
% clc;
% l = 1;
% vol_v = {[ 3*114  2*114  114 32 16 12 8]};
% for vol_i = 1:length(vol_v)
%     [minCurrErr, currSolutionTiling] = loopTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v{vol_i});
% 
% %     [minCurrErr, tilingCurrRes, currSolutionTiling] = recursiveTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v(vol_i), ind2data, tiling, solutionTiling, Inf);
%     if minCurrErr < inf
%         volumeRes(l) = vol_v(vol_i);
%         minErr(l) = minCurrErr;
%         solutionTilingRes(l) = currSolutionTiling;
%         l = l + 1;
%     end
% end
% 
% 
% 
% 
