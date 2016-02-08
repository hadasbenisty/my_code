close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));

overwrite = true;
files = {'tilingDemoData' };
figspath1 = fullfile('tilingDemoOutput');
wrkspname = fullfile(figspath1, 'wrkspace.mat');
if exist(wrkspname, 'file') && ~overwrite
    load(wrkspname);
else
    
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
    % params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trial);
    %% Data
    t = 0:1/1e3:2;
    x=chirp(t, 500, 1, 400, 'quadratic', [] ,'convex');
    P = spectrogram(x,20,4,256,1e3,'yaxis');
    data = log(abs(P(1:7:end,1:7:end)));
    data=data(1:10,1:10);
    figure;
    imagesc(data);colorbar;
    title('Data');
    
    %% Run Qu. 2D
    [row_tree, col_tree] = RunQuestionnaire(params, data);
    
    %% Visualization
    mkNewFolder(figspath1);
    save(wrkspname);
end
% row_aff = CalcEmdAff(data.', col_tree, params.row_emd, params.init_aff_row.on_rows);
% [row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
% [ row_order ] = OrganizeDiffusion3DbyOneDim( data, row_vecs*row_vals );
% [~, row_order1] = sort(row_order);
[~, row_order] = sort(row_tree{2}.clustering);
figure;
plotTreeWithColors(row_tree, (row_order));

[~, col_order] = sort(col_tree{2}.clustering);
figure;
plotTreeWithColors(col_tree, (col_order));

orderedData = data(row_order, :);
% col_aff = CalcEmdAff(data, row_tree, params.col_emd, params.init_aff_row.on_rows);
% [col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
% [ col_order ] = OrganizeDiffusion3DbyOneDim( data, col_vecs*col_vals );
orderedData = data(row_order, :);
orderedData = orderedData(:, col_order);

% orderedData = orderedData(:, col_order);
% prepare trees for recursion
for treeLevel = 1:length(row_tree)
    row_orderedtree{treeLevel} = row_tree{treeLevel};
    row_orderedtree{treeLevel}.clustering = row_tree{treeLevel}.clustering( row_order);
end
for treeLevel = 1:length(col_tree)
    col_orderedtree{treeLevel} = col_tree{treeLevel};
    col_orderedtree{treeLevel}.clustering = col_tree{treeLevel}.clustering( col_order);
end
figure;subplot(2,1,1);
plotTreeWithColors(row_orderedtree, 1:size(data,1));
title('Row Tree');

subplot(2,1,2);
plotTreeWithColors(col_orderedtree, 1:size(data,2))
title('Col Tree');


ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(data, 1), size(data, 2), 1);
tiling.isLeader = zeros(size(data, 1), size(data, 2), 1);
figure;
clc;
for volume = numel(data):-1:2
    minErr = recursiveTiling2D(data(:,:,1), row_orderedtree, col_orderedtree, volume, ind2data, tiling, Inf);
end


% evaluate all trims having constant volume



