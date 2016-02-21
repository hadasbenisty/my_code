close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../gen_utils'));
addpath(genpath('../tiling'));

files = {'8_6_14_1-20_control' '8_6_14_21-60_cno'};
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
dataFor3DQ = X;

%% Run Qu. 3D
params  = SetGenericQuestParamsD8;
% making sure that the timing tree is not too complicated
params.row_tree.treeDepth=5;
[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, permute(dataFor3DQ, [2 3 1]) );
figure;
subplot(1,2,1);
[vecs_t, vals] = CalcEigs(threshold(col_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs_t * vals, 1:size(dataFor3DQ, 2), 'Time Embedding');
subplot(1,2,2);
plotEmbeddingWithColors(vecs_t * vals, [ones(1, 100) 100*ones(1, 20) 200*ones(1, 240)], 'Time Colored by Tone');

figure;
[vecs_r, vals] = CalcEigs(threshold(row_dual_aff, 0.0), 3);
plotEmbeddingWithColors(vecs_r * vals, 1:size(dataFor3DQ, 1), 'Nuerons Embedding');
figure;
[vecs_T, vals] = CalcEigs(threshold(trial_dual_aff, 0.0), 4);
plotEmbeddingWithColors(vecs_T * vals, 1:size(dataFor3DQ, 3), 'Trials Embedding');

figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trial_tree, 1:length(trial_dual_aff));    title('Trialsl Tree');

%% Order by tree
meanData = mean(dataFor3DQ, 3);
[~, row_order] = sort(row_tree{2}.clustering);
[~, col_order] = sort(col_tree{2}.clustering);
orderedData = meanData(row_order, :);
orderedData = orderedData(:, col_order);

% prepare trees for recursion
for treeLevel = 1:length(row_tree)
    row_orderedtree{treeLevel} = row_tree{treeLevel};
    row_orderedtree{treeLevel}.clustering = row_tree{treeLevel}.clustering( row_order);
end
for treeLevel = 1:length(col_tree)
    col_orderedtree{treeLevel} = col_tree{treeLevel};
    col_orderedtree{treeLevel}.clustering = col_tree{treeLevel}.clustering( col_order);
end

figure;
subplot(2,2,1);
imagesc(meanData);title('Orig Data');
xlabel('Time');
set(gca, 'YTick', 1:length(NeuronsLabels));
set(gca, 'YTickLabel',NeuronsLabels);
subplot(2,2,4);
imagesc(orderedData);title('Ordered Data');
xlabel('Time');
set(gca, 'YTick', 1:length(NeuronsLabels));
set(gca, 'YTickLabel',NeuronsLabels(row_order));

subplot(2,2,3);
plotTreeWithColors(row_orderedtree, 1:size(meanData,1));
title('Row Tree'); colorbar('off')
view(-90,90)
subplot(2,2,2);
plotTreeWithColors(col_orderedtree, 1:size(meanData,2))
title('Col Tree');colorbar('off')


ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
solutionTiling = [];
figure;
clc;
l = 1;
vol_v = {[  1e3:-1:361 359:-1:109 107:-1:100  ]};
for vol_i = 1:length(vol_v)
    [minCurrErr, currSolutionTiling] = loopTiling2D(orderedData, row_orderedtree, col_orderedtree, vol_v{vol_i});

%     [minCurrErr, tilingCurrRes, currSolutionTiling] = recursiveTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v(vol_i), ind2data, tiling, solutionTiling, Inf);
    if minCurrErr < inf
        volumeRes(l) = vol_v(vol_i);
        minErr(l) = minCurrErr;
        solutionTilingRes(l) = currSolutionTiling;
        l = l + 1;
    end
end

[~,a]=sort(col_order);
for ci = 1:length(unique(solutionTilingRes.isbusy(:)))
   [ind_i, ind_j] = find(solutionTilingRes.isbusy == ci) ;
   meanTiled(ind_i, ind_j) = mean(mean(orderedData(ind_i, ind_j)));
end

%% plot data and tiling of the other experiments

X1 = X(:, :, 41:75);
X2 = X(:, :, 76:120);
X3 = X(:, :, 121:end);
meanData1 = mean(X1, 3);
meanData2 = mean(X2, 3);
meanData3 = mean(X3, 3);

orderedData1 = meanData1(row_order, :);
orderedData1 = orderedData1(:, col_order);
orderedData2 = meanData2(row_order, :);
orderedData2 = orderedData2(:, col_order);
orderedData3 = meanData3(row_order, :);
orderedData3 = orderedData3(:, col_order);
for ci = 1:length(unique(solutionTilingRes.isbusy(:)))
   [ind_i, ind_j] = find(solutionTilingRes.isbusy == ci) ;
   meanTiled1(ind_i, ind_j) = mean(mean(orderedData1(ind_i, ind_j)));
   meanTiled2(ind_i, ind_j) = mean(mean(orderedData2(ind_i, ind_j)));
   meanTiled3(ind_i, ind_j) = mean(mean(orderedData3(ind_i, ind_j)));
end
figure;plotTiledData(orderedData(:,a), meanTiled(:,a), NeuronsLabels(row_order), solutionTilingRes.isbusy(:,a), '8/12/14')
figure;plotTiledData(orderedData1(:,a), meanTiled1(:,a), NeuronsLabels(row_order), solutionTilingRes.isbusy(:,a), '8/15/13')
figure;plotTiledData(orderedData2(:,a), meanTiled2(:,a), NeuronsLabels(row_order), solutionTilingRes.isbusy(:,a), '8/17/14 part 1')
figure;plotTiledData(orderedData3(:,a), meanTiled3(:,a), NeuronsLabels(row_order), solutionTilingRes.isbusy(:,a), '8/17/14 part 2')




