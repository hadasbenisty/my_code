close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../Questionnaire'));
addpath(genpath('../utils'));

rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
eigsnum_col = 12;
eigsnum_row = 12;
eigsnum_trial= 12;
row_alpha = 0;
row_beta = 0;
col_alpha = 0;
col_beta = 0;
trial_alpha = 0;
trial_beta = 0;
% params  = SetQuest3DParams(eigsnum_col, eigsnum_row, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trial);
%% Load Data

datapth = '..\..\..\datasets\biomed\D8';

files = {'8_6_14_1-20_control' '8_6_14_21-60_cno'};
nt = 120;
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);


[nr, nt, nT] = size(X);
% perm data

data = X;
%% Run Qu. 3D
[row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data);


figure;
subplot(3,2,1);
[vecs, vals] = CalcEigs(threshold(col_aff, 0.2), 2);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(3,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1, 40) 100*ones(1, 80)], 'Time Colored by Tone');

subplot(3,1,2);
[vecs, vals] = CalcEigs(threshold(row_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
subplot(3,2,5);
[vecs, vals] = CalcEigs(threshold(trials_aff, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');
subplot(3,2,6);
plotEmbeddingWithColors(vecs * vals, [ones(1, 20) 100*ones(1, 40)], 'Trials Colored by CNO');
figure;    subplot(3,1,1);
plotTreeWithColors(col_tree, 1:length(col_aff));    title('Time Col Tree');
subplot(3,1,2);    plotTreeWithColors(row_tree, 1:length(row_aff));    title('Nuerons Tree')
subplot(3,1,3);    plotTreeWithColors(trials_tree, 1:length(trials_aff));    title('Trialsl Tree');

tree_level=4;
% clusters = unique(row_tree{tree_level}.clustering);
% for ci = 1:length(clusters)
%    inds2neurons = (find(row_tree{tree_level}.clustering == clusters(ci)));
%    N = length(inds2neurons);
%    R = ceil(sqrt(N));
%    figure;
%    for r = 1:N
%       subplot(R, R, r);
%       imagesc(permute(data(inds2neurons(r), :, :), [3 2 1]));      
%    end
%     
% end

[vecs, vals] = CalcEigs(threshold(col_aff, 0.4), 2);
[~,toneIndic] = min(diff(vecs(2:end)));
[vecs, vals] = CalcEigs(threshold(trials_aff, 0.0), 2);
CNO_th = min(vecs(1:20));
inds2CNO = find((vecs )<= CNO_th);
inds2CONTROL = find(vecs > CNO_th);
inds2CONTROL_part1 = inds2CONTROL(inds2CONTROL<=30);
inds2CONTROL_part2 = inds2CONTROL(inds2CONTROL>30);
[meanMat_beforeCNO_beforeTone, allMat_beforeCNO_beforeTone] = getCentroidsByTree(row_tree{tree_level}, X(:,1:toneIndic,inds2CONTROL_part1), NeuronsLabels, NeuronsLabels);
[meanMat_afterCNO_beforeTone, allMat_afterCNO_beforeTone] = getCentroidsByTree(row_tree{tree_level}, X(:,1:toneIndic,inds2CONTROL_part2), NeuronsLabels, NeuronsLabels);
[meanMat_CNO_beforeTone, allMat_CNO_beforeTone] = getCentroidsByTree(row_tree{tree_level}, X(:,1:toneIndic,inds2CNO), NeuronsLabels, NeuronsLabels);

[ aff_mat1,  ] = CalcInitAff( meanMat_beforeCNO_beforeTone.', params.init_aff_row );
[vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);% 
[ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat_beforeCNO_beforeTone, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order1 ,:) * vals, 1:13, 'Nuerons Embedding');
plotByClustering(meanMat_beforeCNO_beforeTone(row_order1, :), allMat_beforeCNO_beforeTone(row_order1), ['Before Tone Before CNO Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat_afterCNO_beforeTone(row_order1, :), allMat_afterCNO_beforeTone(row_order1), ['Before Tone After CNO Organized By Tree Of Before CNO Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat_CNO_beforeTone(row_order1, :), allMat_CNO_beforeTone(row_order1), ['Before Tone During CNO Organized By Tree Of Before CNO Tree Level ' num2str(tree_level)]);




[meanMat_beforeCNO_afterTone, allMat_beforeCNO_afterTone] = getCentroidsByTree(row_tree{tree_level}, X(:,toneIndic+1:end,inds2CONTROL_part1), NeuronsLabels, NeuronsLabels);
[meanMat_afterCNO_afterTone, allMat_afterCNO_afterTone] = getCentroidsByTree(row_tree{tree_level}, X(:,toneIndic+1:end,inds2CONTROL_part2), NeuronsLabels, NeuronsLabels);
[meanMat_CNO_afterTone, allMat_CNO_afterTone] = getCentroidsByTree(row_tree{tree_level}, X(:,toneIndic+1:end,inds2CNO), NeuronsLabels, NeuronsLabels);

% [ aff_mat1,  ] = CalcInitAff( meanMat_beforeCNO_afterTone.', params.init_aff_row );
% [vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);% 
% [ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat_before, vecs*vals );
% figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
plotByClustering(meanMat_beforeCNO_afterTone(row_order1, :), allMat_beforeCNO_afterTone(row_order1), ['After Tone Before CNO Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat_afterCNO_afterTone(row_order1, :), allMat_afterCNO_afterTone(row_order1), ['After Tone After CNO Organized By Tree Of Before CNO Tree Level ' num2str(tree_level)]);
plotByClustering(meanMat_CNO_afterTone(row_order1, :), allMat_CNO_afterTone(row_order1), ['After Tone During CNO Organized By Tree Of Before CNO Tree Level ' num2str(tree_level)]);


figure;
imagesc([meanMat_beforeCNO_beforeTone(row_order1, :) meanMat_beforeCNO_afterTone(row_order1, :);...
         meanMat_CNO_beforeTone(row_order1, :) meanMat_CNO_afterTone(row_order1, :);...
         meanMat_afterCNO_beforeTone(row_order1, :) meanMat_afterCNO_afterTone(row_order1, :)])
[X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_4_14_1-25_control'}, nt);
[meanMat1_beforeTone, allMat1_beforeTone] = getCentroidsByTree(row_tree{tree_level}, X1(:,1:toneIndic,:), NeuronsLabels, NeuronsLabels1);
[meanMat1_afterTone, allMat1_afterTone] = getCentroidsByTree(row_tree{tree_level}, X1(:,toneIndic+1:end,:), NeuronsLabels, NeuronsLabels1);
figure;
imagesc([meanMat_beforeCNO_beforeTone(row_order1, :) meanMat_beforeCNO_afterTone(row_order1, :);...
    meanMat1_beforeTone(row_order1, :) meanMat1_afterTone(row_order1, :)]);

%   [X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'7_23_14_1-35_control'}, nt);

    
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