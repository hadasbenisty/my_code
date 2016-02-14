close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));

files = { '8_12_14_1-40'  '8_15_13_1-35' '8_17_14_1-45'  '8_17_14_46-80'};% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
eigsnum_col = 12;eigsnum_row = 12;eigsnum_trial= 12;
row_alpha = .2;row_beta = 0;col_alpha = .2;col_beta = 0;trial_alpha = .2;trial_beta = 0;
params  = SetQuest3DParams(eigsnum_row, eigsnum_col, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
%% Load Data

datapth = '..\..\..\datasets\biomed\D30';

nt = 360;
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);


[nr, nt, nT] = size(X(:, :, 1:40));
data = X(:, :, 1:40);

%% Run Qu. 3D
params  = SetGenericQuestParams;
   
[ row_tree, col_tree, trial_tree, row_dual_aff, col_dual_aff, trial_dual_aff ] = RunGenericQuestionnaire3D( params, data );
figure;
subplot(3,2,1);
[vecs, vals] = CalcEigs(threshold(col_dual_aff, 0.2), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(3,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 100) 100*ones(1, 20) 200*ones(1, 240)], 'Time Colored by Tone');

subplot(3,1,2);
[vecs, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
subplot(3,1,3);
[vecs, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');
% subplot(3,2,6);
% plotEmbeddingWithColors(vecs * vals, [ones(1, 35) 100*ones(1, 40) 150*ones(1, 45) 200*ones(1, 35) ], 'Trials Colored By Experiment');

figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trial_tree, 1:length(trial_dual_aff));    title('Trialsl Tree');

meanMat=[];l = 2;
folders = unique(row_tree{l}.clustering);
for ci = 1:length(folders)
    inds2folder = find(row_tree{l}.clustering == folders(ci));
    meanMat(ci, :) = mean(mean(permute(data(inds2folder, :, :), [2 3 1]),2), 3);
end
figure;
imagesc(meanMat);
set(gca, 'Ytick', 1:length(folders));


X1 = X(:, :, 41:75);
X2 = X(:, :, 76:120);
X3 = X(:, :, 121:end);
meanMat=[];l = 2;
folders = unique(row_tree{l}.clustering);
for ci = 1:length(folders)
    inds2folder = find(row_tree{l}.clustering == folders(ci));
    meanMat(ci, :) = mean(mean(permute(X1(inds2folder, :, :), [2 3 1]),2), 3);
end
figure;
imagesc(meanMat);
title('8_15_13_1-35');
set(gca, 'Ytick', 1:length(folders));
meanMat=[];l = 2;
folders = unique(row_tree{l}.clustering);
for ci = 1:length(folders)
    inds2folder = find(row_tree{l}.clustering == folders(ci));
    meanMat(ci, :) = mean(mean(permute(X2(inds2folder, :, :), [2 3 1]),2), 3);
end
figure;
imagesc(meanMat);
title('8_17_14_1-45');
set(gca, 'Ytick', 1:length(folders));
meanMat=[];l = 2;
folders = unique(row_tree{l}.clustering);
for ci = 1:length(folders)
    inds2folder = find(row_tree{l}.clustering == folders(ci));
    meanMat(ci, :) = mean(mean(permute(X3(inds2folder, :, :), [2 3 1]),2), 3);
end
figure;
imagesc(meanMat);
title('8_17_14_46-80');
set(gca, 'Ytick', 1:length(folders));
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
