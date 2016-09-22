close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));
addpath(genpath('../../grangerCausality/'));
experiment = 'M2';
dims4quest = 3;

datapth = '..\..\..\datasets\biomed\M2';
nt = 360;
ind2data =3;
toneVec = [ones(1,40) 100*ones(1, 80)];
files = { '4_4_14' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'


rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
%% Load Data


[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);




data = X(:, ind2data:3:end, :);
% selectedTimeFrams = 35:60;

selectedTimeFrams=1:size(data, 2);
data = data(:, selectedTimeFrams, :);
data_part1 = data(1:47, :, :);
data_part2 = data(48:end, :, :);
[nr, nt, nT] = size(data);
toneVec=toneVec(selectedTimeFrams);
neuron_tree_level = 2;
istopdown = true;

%% Run Qu.
params = SetGenericDimsQuestParams(dims4quest, istopdown);
if istopdown
    params.data.over_rows = false;
    runningOrder = [1 2 3];
    for ind = 1:dims4quest
        params.emd{ind}.beta = .1;
        params.tree{ind}.runOnEmbdding = false;
        params.tree{ind}.splitsNum = 2;
        params.tree{ind}.treeDepth = inf;
        
    end
params.tree{1}.runOnEmbdding = true;
    params.init_aff{3}.thresh = .1;
    
else
    params.data.over_rows = false;
    runningOrder = [2 3 1];
    params.init_aff{1}.thresh = 0.0;
    params.init_aff{2}.thresh = 0;
    params.init_aff{3}.thresh = 0.0;
    params.emd{1}.beta = 0.1;
    params.emd{2}.beta = 0.2;
    params.emd{3}.beta = 0.2;

end


% Run Qu.
% params = SetGenericDimsQuestParams(dims4quest, true);
% for ind = 1:dims4quest
%     params.emd{ind}.beta = 0.5;
%     params.tree{ind}.splitsNum = 2;
%     params.tree{ind}.treeDepth = inf;
%     params.tree{ind}.runOnEmbdding = false;
% end
% params.tree{1}.runOnEmbdding = true;
% params.tree{2}.runOnEmbdding = false;
% params.tree{3}.runOnEmbdding = false;
% if dims4quest == 3
%     runningOrder = [1 2 3];
% else
%     runningOrder = [1  2];
% end
params.data.to_normalize = true;
params.data.normalization_type = 'by_std';
% params.init_aff{3}.metric = 'euc';
% params.init_aff{3}.metric = 'cosine_similarityOnTrials';
params.init_aff{3}.metric = 'cosine_similarity';
params.n_iters = 2;

params.verbose = 2;

[ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params, permute(data,(runningOrder) ) );
[ Trees_part1, dual_aff_part1, init_aff_part1 ] = RunGenericDimsQuestionnaire( params, permute(data_part1,(runningOrder) ) );
[ Trees_part2, dual_aff_part2, init_aff_part2 ] = RunGenericDimsQuestionnaire( params, permute(data_part2,(runningOrder) ) );

% measure the accuracy for neurons level
clusters = unique(Trees{runningOrder==1}{end-1}.clustering);
for ci = 1:length(clusters)
    inds = find(Trees{runningOrder==1}{end-1}.clustering==ci);
    Acc(ci) = numel(intersect(inds, 1:46))/numel(inds);
    
end
Acc = max(Acc, 1-Acc);
figure;plotTreeWithColors(Trees{runningOrder==1}, [ones(1,46) 100*ones(1, 35)]);
title(['Neurons Tree, Colored By Layer Acc. = ' num2str(Acc)]);


figure;
[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==2}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time - ');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Time Colored By Tone - ');
figure;subplot(2,2,1);
[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==1}, 0.25), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Neurons - ');
subplot(2,2,2);
plotEmbeddingWithColors(vecs * vals, [ones(1,46) 100*ones(1, 35)], 'Nuerons Colored by Layer - ');
subplot(2,2,3);
[vecs, vals] = CalcEigs(threshold(dual_aff_part1{runningOrder==1}, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data_part1, 1), 'Part 1: Neurons - ');
subplot(2,2,4);
[vecs, vals] = CalcEigs(threshold(dual_aff_part2{runningOrder==1}, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data_part2, 1), 'Part 2: Neurons - ');


figure;
[vecs, vals] = CalcEigs(threshold(dual_aff_part1{runningOrder==2}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Part 1: Time - ');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Part 1: Time Colored By Tone - ');

figure;
[vecs, vals] = CalcEigs(threshold(dual_aff_part2{runningOrder==2}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Part 2: Time - ');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Part 2: Time Colored By Tone - ');


if dims4quest==3
    figure; subplot(1, 3, 1);
    [vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==3}, 0.0), 3);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials - ');
    subplot(1, 3, 2);
    [vecs, vals] = CalcEigs(threshold(dual_aff_part1{runningOrder==3}, 0.0), 3);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Part 1: Trials - ');
    subplot(1, 3, 3);
    [vecs, vals] = CalcEigs(threshold(dual_aff_part2{runningOrder==3}, 0.0), 3);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Part 2: Trials - ');
end

%% ordering by nuerons clusters - A
% selectedTimeFrams = 1:120;
dataN = NormalizeData(data, params.data);
dataN_part1 = NormalizeData(data_part1, params.data);
dataN_part2 = NormalizeData(data_part2, params.data);

[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(Trees{runningOrder==1}{neuron_tree_level}, dataN, NeuronsLabels, NeuronsLabels);
affine = feval(params.init_aff{runningOrder==1}.initAffineFun, permute(meanMatAlltrials, [2 1 3]), params.init_aff{1});
[vecs, vals] = CalcEigs(affine, 2);
classes_order = OrganizeDiffusion3DbyOneDim( meanMatAlltrials, vecs*vals );
normalized_centroids = plotByClustering(allMat(classes_order),  'All Data', selectedTimeFrams);%selectedTimeFrams
   
[meanMat_part1, allMat_part1, meanMatAlltrials_part1] = getCentroidsByTree(Trees_part1{runningOrder==1}{neuron_tree_level}, dataN_part1, NeuronsLabels, NeuronsLabels);
affine_part1 = feval(params.init_aff{runningOrder==1}.initAffineFun, permute(meanMatAlltrials_part1, [2 1 3]), params.init_aff{1});
[vecs, vals] = CalcEigs(affine_part1, 2);
classes_order_part1 = OrganizeDiffusion3DbyOneDim( meanMatAlltrials_part1, vecs*vals );
normalized_centroids_part1 = plotByClustering(allMat_part1(classes_order_part1),  'Part 1', selectedTimeFrams);%selectedTimeFrams
 
[meanMat_part2, allMat_part2, meanMatAlltrials_part2] = getCentroidsByTree(Trees_part2{runningOrder==1}{neuron_tree_level}, dataN_part2, NeuronsLabels, NeuronsLabels);
affine_part2 = feval(params.init_aff{runningOrder==1}.initAffineFun, permute(meanMatAlltrials_part2, [2 1 3]), params.init_aff{1});
[vecs, vals] = CalcEigs(affine_part2, 2);
classes_order_part2 = OrganizeDiffusion3DbyOneDim( meanMatAlltrials_part2, vecs*vals );
normalized_centroids_part2 = plotByClustering(allMat_part2(classes_order_part2),  'Part 2', selectedTimeFrams);%selectedTimeFrams

plotOnsetByClustering(normalized_centroids,  'All Data', 1:size(normalized_centroids, 2), .8,false);
plotOnsetByClustering(normalized_centroids_part1,  'Part 1', 1:size(normalized_centroids, 2), .8,false);
plotOnsetByClustering(normalized_centroids_part2,  'Part 2', 1:size(normalized_centroids, 2), .8,false);

p.alpha = .7;
p.verbose = 2;
p.momax = 10;% maximum model order for model order estimation
p.icregmode = 'OLS';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
p.morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
p.regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
p.tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
p.mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
p.acmaxlags = 1500;
p.useSS = true;
p.alpha = .005;    [F, pval, sig] = grangerCausalityWrapper(permute(meanMatAlltrials(classes_order, :, :), [1 2 3]), p);
p.alpha = .35;    [F, pval, sig] = grangerCausalityWrapper([meanMatAlltrials_part1(classes_order_part1, :, :); meanMatAlltrials_part2(classes_order_part2, :, :)], p);


% [F1,c_v1] = granger_causeWrapper(meanMat(classes_order,:),0.0005,40);
% [F1,c_v1] = granger_causeWrapper([meanMat_part1(classes_order_part1, :); meanMat_part2(classes_order_part2, :)],0.00001,40);

% figure;imagesc(F1<=c_v1);colormap gray
% getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, :, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' By Tree Level ' num2str(neuron_tree_level)])
%        
% getClusteringByTreeAndPlot(Trees_part1{runningOrder==1}{neuron_tree_level}, data_part1(:, :, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' - Part 1 - By Tree Level ' num2str(neuron_tree_level)])
% 
% getClusteringByTreeAndPlot(Trees_part2{runningOrder==1}{neuron_tree_level}, data_part2(:, :, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' - Part 2 - By Tree Level ' num2str(neuron_tree_level)])

% selectedTimeFrams = 35:60;
% getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' By Tree Level ' num2str(neuron_tree_level)])
%        
% getClusteringByTreeAndPlot(Trees_part1{runningOrder==1}{neuron_tree_level}, data_part1(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' - Part 1 - By Tree Level ' num2str(neuron_tree_level)])
% 
% getClusteringByTreeAndPlot(Trees_part2{runningOrder==1}{neuron_tree_level}, data_part2(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} ,...
%     [experiment ' - Part 2 - By Tree Level ' num2str(neuron_tree_level)])




% 
% 
% %% ordering neurons clusters - B
% clusteringN = Trees{runningOrder==1}{neuron_tree_level}.clustering;
% [~, Norder] = sort(clusteringN);
% clusteringt = Trees{runningOrder==2}{2}.clustering;
% clusteringT = Trees{runningOrder==3}{2}.clustering;
% 
% [~, torder] = sort(clusteringt);
% [~, Torder] = sort(clusteringT);
% clstrsN = unique(clusteringN);
% clstrst = unique(clusteringt);
% 
% clusteredData = data(Norder, :, :);
% 
% figure;imagesc(mean(clusteredData(:, :,:), 3));
% y = -.50;x = -.50;
% for r = 1:length(clstrsN)
%     line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     
%     line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
% end
% title([experiment ' - Ordering by Neurons '  'Averaged Over Trials ']);
% 
% 
% %% %% Part 1 ordering neurons clusters
% clusteringN = Trees_part1{runningOrder==1}{neuron_tree_level}.clustering;
% [~, Norder] = sort(clusteringN);
% clusteringt = Trees_part1{runningOrder==2}{2}.clustering;
% clusteringT = Trees_part1{runningOrder==3}{2}.clustering;
% 
% [~, torder] = sort(clusteringt);
% [~, Torder] = sort(clusteringT);
% clstrsN = unique(clusteringN);
% clstrst = unique(clusteringt);
% 
% clusteredData_part1 = data_part1(Norder, :, :);
% 
% figure;imagesc(mean(clusteredData_part1(:, :,:), 3));
% y = -.50;x = -.50;
% for r = 1:length(clstrsN)
%     line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     
%     line([1 size(data_part1, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     y = y + Trees_part1{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
% end
% title([experiment ' - Part 1 - Ordering by Neurons '  'Averaged Over Trials ']);
% 
% %% %% Part 2 ordering neurons clusters
% clusteringN = Trees_part2{runningOrder==1}{neuron_tree_level}.clustering;
% [~, Norder] = sort(clusteringN);
% clusteringt = Trees_part2{runningOrder==2}{2}.clustering;
% clusteringT = Trees_part2{runningOrder==3}{2}.clustering;
% 
% [~, torder] = sort(clusteringt);
% [~, Torder] = sort(clusteringT);
% clstrsN = unique(clusteringN);
% clstrst = unique(clusteringt);
% 
% clusteredData_part2 = data_part2(Norder, :, :);
% 
% figure;imagesc(mean(clusteredData_part2(:, :,:), 3));
% y = -.50;x = -.50;
% for r = 1:length(clstrsN)
%     line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     
%     line([1 size(data_part2, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%     y = y + Trees_part2{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
% end
% title([experiment ' - Part 2 - Ordering by Neurons '  'Averaged Over Trials ']);
