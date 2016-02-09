close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));


figspath1 = fullfile('tilingDemoOutput');


rng(654164);

%% Init params

eigsnum_col = 2;
eigsnum_row = 3;
eigsnum_trial= 2;
row_alpha = .2;
row_beta = 0;
col_alpha = .2;
col_beta = 0;
trial_alpha = .2;
trial_beta = 0;
params  = SetQuest3DParams(eigsnum_col, eigsnum_row, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
params.init_aff_row_metric = 'euc';
params.init_aff_col_metric = 'euc';
%% Data

data = [1 1  1   1  2  2  3  3  3  3;...
    1 1  1   1  2  2  3  3  3  3;...
    4 4  8   8  11 11 11 11 80 80;...
    4 4  8   8  12 12 12  12 80 80;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 70 70 70 70 90 90;...
    6 14 16  17 95 95 95 95 90 90];
data = data + randn(size(data))*0.3;
figure;
imagesc(data);colorbar;
title('Data');

%% Run Qu. 2D
[row_tree, col_tree] = RunQuestionnaire(params, data);

%% Visualization

[~, row_order] = sort(row_tree{2}.clustering);
[~, col_order] = sort(col_tree{2}.clustering);
orderedData = data(row_order, :);
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
imagesc(data);title('Data');
subplot(2,2,2);
imagesc(orderedData);title('Ordered Data');
subplot(2,2,3);
plotTreeWithColors(row_orderedtree, 1:size(data,1));
title('Row Tree');

subplot(2,2,4);
plotTreeWithColors(col_orderedtree, 1:size(data,2))
title('Col Tree');


ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
solutionTiling = [];
figure;
clc;
l = 1;
vol_v = numel(orderedData):-1:4;

for vol_i = 1:length(vol_v)
    [minCurrErr, tilingCurrRes, currSolutionTiling] = recursiveTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v(vol_i), ind2data, tiling, solutionTiling, Inf);
    [minCurrErr1, currSolutionTiling1] = loopTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v(vol_i));
    
    if minCurrErr1 ~= minCurrErr
        if minCurrErr~= inf
            if max(max(abs(currSolutionTiling.isbusy - currSolutionTiling1.isbusy)))~=0
                error('s');
            end
        end
    end
    if minCurrErr < inf
        volumeRes(l) = vol_v(vol_i);
        minErr(l) = minCurrErr;
        tilingRes(l) =   tilingCurrRes;
        solutionTilingRes(l) = currSolutionTiling;
        l = l + 1;
    end
end
volumeRes(end+1) = 1;
minErr(end+1) = 0;
figure;plot(volumeRes, minErr, '-*');

%
%
% small_row_tree{1}.folder_count = 4;
% small_row_tree{1}.folder_sizes = [1 1 1 1];
% small_row_tree{1}.clustering = [1 2 3 4];
% small_row_tree{1}.super_folders = [1 1 2 2];
%
% small_row_tree{2}.folder_count = 2;
% small_row_tree{2}.folder_sizes = [2 2];
% small_row_tree{2}.clustering = [1 1 2 2];
% small_row_tree{2}.super_folders = [ 1 1];
%
% small_row_tree{3}.folder_count = 1;
% small_row_tree{3}.folder_sizes = 4;
% small_row_tree{3}.clustering = [1 1 1 1];
% small_row_tree{3}.super_folders = [];
%
% small_col_tree = small_row_tree;
% small_data = rand(4);
% small_vol = 4;
% ind2data = 1;
% minErr = Inf;
% tiling.isbusy = zeros(size(small_data, 1), size(small_data, 2), 1);
% tiling.isLeader = zeros(size(small_data, 1), size(small_data, 2), 1);
% solutionTiling = [];
% figure;
% clc;
% [minCurrErr, tilingCurrRes, currSolutionTiling] = recursiveTiling2D(small_data, small_row_tree, small_col_tree, small_vol, 1, tiling, solutionTiling, Inf);
% [minErr1, solutionTiling] = loopTiling2D(small_data, small_row_tree, small_col_tree, small_vol);
