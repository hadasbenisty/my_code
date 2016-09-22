close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));

addpath(genpath('../tiling'));


figspath1 = fullfile('tilingDemoOutput');


rng(654164);

%% Data

X = [1 1  1   1  2  2  3  3  3  3;...
    1 1  1   1  2  2  3  3  3  3;...
    4 4  8   8  11 11 11 11 80 80;...
    4 4  8   8  12 12 12  12 80 80;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 70 70 70 70 90 90;...
    6 14 16  17 95 95 95 95 90 90];
data = X + randn(size(X))*0.3;
figure;
imagesc(data);colorbar;
title('Data');

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
% params  = SetQuestPCAclusteringParams;
params = SetGenericDimsQuestParams(2, true);
params.tree{1}.splitsNum=5;
params.tree{2}.splitsNum=5;
% params.tree{1}.clusteringAlgo = @MTSGClassWrapper;
% params.tree{1}.min_cluster = 1;
% params.tree{2}.clusteringAlgo = @MTSGClassWrapper;
% params.tree{2}.min_cluster = 1;
[ Trees, dual_aff ] = RunGenericDimsQuestionnaire( params, data  );
%% Visualization

[~, row_order] = sort(Trees{1}{2}.clustering);
[~, col_order] = sort(Trees{2}{2}.clustering);
orderedData = data(row_order, :);
orderedData = orderedData(:, col_order);

% prepare trees for recursion
for treeLevel = 1:length(Trees{1})
    row_orderedtree{treeLevel} = Trees{1}{treeLevel};
    row_orderedtree{treeLevel}.clustering = Trees{1}{treeLevel}.clustering( row_order);
end
for treeLevel = 1:length(Trees{2})
    col_orderedtree{treeLevel} = Trees{2}{treeLevel};
    col_orderedtree{treeLevel}.clustering = Trees{2}{treeLevel}.clustering( col_order);
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
classes = unique(X);
for ci=1:length(classes)
    tilesSizes(ci) = length(find(X(:)==classes(ci)));
end
p1.vol = unique(tilesSizes);
p1.vol = sort(p1.vol, 'ascend');

p1.beta_col = 0;
p1.beta_row = 0;
p1.verbose = 2;
%             p1.err_fun = @evalTilingErr;
p1.err_fun = @evalTilingErrTauMeas;p1.tau = 0.1;
tiling.isbusy=zeros(size(data));


[minCurrErr, currSolutionTiling] = loopTiling2DEfficient(orderedData, row_orderedtree, col_orderedtree, p1);
[minCurrErr1, tilingCurrRes1, currSolutionTiling1] = recursiveTiling2D(orderedData, row_orderedtree, col_orderedtree, p1.vol, 1, tiling, [], inf);

[~, meanTiled] = evalTilingErr(X, currSolutionTiling, p1);
figure;plotTiledData(orderedData, meanTiled, '', currSolutionTiling.isbusy, [ ' \tau = ' num2str(p1.tau)])

[~, meanTiled1] = evalTilingErr(X, currSolutionTiling1, p1);
figure;plotTiledData(orderedData, meanTiled1, '', currSolutionTiling1.isbusy, [ ' \tau = ' num2str(p1.tau)])

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
