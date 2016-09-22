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

% files = { '8_12_14_1-40' '8_15_13_1-35'  '8_17_14_1-45'  '8_17_14_46-80'};% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'
datapth = '..\..\..\datasets\biomed\D30';
nt = 360;
toneVec = [ones(1,120) 100*ones(1, 240)];
files = { '8_12_14_1-40' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

nt = 360;
toneVec = [ones(1,120) 100*ones(1, 240)];
files = { '8_12_14_1-40' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'


%% Load Data


[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);

params.tree{1}.runOnEmbdding = true;
params.tree{2}.runOnEmbdding = true;%false;
params.tree{3}.runOnEmbdding = true;%false;
selectedTimeFrams = 115:140;
data_all = X;
data_tone =  X(:, selectedTimeFrams, :);

%% Run Qu.
dims4quest=3;
params = SetGenericDimsQuestParams(dims4quest, true);
for ind = 1:dims4quest
    params.emd{ind}.beta = .1;  
    params.init_aff{ind}.metric = 'cosine_similarity';

end

if dims4quest == 3
    runningOrder = [  1 2 3];
else
    runningOrder = [1  2];
end
params.data.over_rows = true;
params.data.to_normalize = true;%false
params.data.normalization_type = 'by_std';

params.n_iters = 1;
params.tree{1}.runOnEmbdding = false;
params.tree{2}.runOnEmbdding = false;
params.tree{3}.runOnEmbdding = false;
params.tree{1}.splitsNum = 10;
params.tree{1}.treeDepth = 4;

params.tree{2}.splitsNum = 2;
params.tree{2}.treeDepth = 6;

params.tree{3}.splitsNum = 5;
params.tree{3}.treeDepth = 3;
params.verbose = 2;





[ Trees_all, dual_aff_all, init_aff_all ] = RunGenericDimsQuestionnaire( params, permute(data_all,(runningOrder) ) );
% [ Trees_tone, dual_aff_tone, init_aff_tone ] = RunGenericDimsQuestionnaire( params, permute(data_tone,(runningOrder) ) );


figure;
[vecs, vals] = CalcEigs(threshold(dual_aff_all{runningOrder==2}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data_all, 2), 'Time - ');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Time Colored By Tone - ');

figure;
[vecs, vals] = CalcEigs(threshold(dual_aff_all{runningOrder==1}, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data_all, 1), 'Nuerons - ');

if dims4quest==3
    figure;
    [vecs, vals] = CalcEigs(threshold(dual_aff_all{runningOrder==3}, 0.0), 2);
    plotEmbeddingWithColors(vecs * vals, 1:size(data_all, 3), 'Trials - ');
end

[X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, 360);
[X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, 360);
[X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, 360);



dataN_all = NormalizeData(data_all, params.data);
dataN_all1 = NormalizeData(X1, params.data);
dataN_all2 = NormalizeData(X2, params.data);
dataN_all3 = NormalizeData(X3, params.data);

dataN_tone = NormalizeData(data_tone, params.data);
neuron_tree_level = 2;
[meanMat_all, allMat_all, meanMatAlltrials_all] = getCentroidsByTree(Trees_all{runningOrder==1}{neuron_tree_level}, dataN_all, NeuronsLabels, NeuronsLabels);
[meanMat_all1, allMat_all1, meanMatAlltrials_all1] = getCentroidsByTree(Trees_all{runningOrder==1}{neuron_tree_level}, dataN_all1, NeuronsLabels, NeuronsLabels1);
[meanMat_all2, allMat_all2, meanMatAlltrials_all2] = getCentroidsByTree(Trees_all{runningOrder==1}{neuron_tree_level}, dataN_all2, NeuronsLabels, NeuronsLabels2);
[meanMat_all3, allMat_all3, meanMatAlltrials_all3] = getCentroidsByTree(Trees_all{runningOrder==1}{neuron_tree_level}, dataN_all3, NeuronsLabels, NeuronsLabels3);
affine_all = feval(params.init_aff{runningOrder==1}.initAffineFun, permute(meanMatAlltrials_all, [2 1 3]), params.init_aff{1});
[vecs, vals] = CalcEigs(affine_all, 10);
classes_order_all = OrganizeDiffusion3DbyOneDim( meanMatAlltrials_all, vecs*vals );
normalized_centroids_all = plotByClustering(allMat_all(classes_order_all),  ['8/12/14 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all1 = plotByClustering(allMat_all1(classes_order_all),  ['8/15/13 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all2 = plotByClustering(allMat_all2(classes_order_all),  ['8/17/14 part 1 Organized By Tree Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all3 = plotByClustering(allMat_all3(classes_order_all),  ['8/17/14 part 2 Organized By Tree Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)]);%selectedTimeFrams

%% Order by tree
meanData = mean(dataN_all, 3);
[~, neurons_order] = sort(Trees_all{runningOrder==1}{2}.clustering);
[~, time_order] = sort(Trees_all{runningOrder==2}{2}.clustering);
orderedData = meanData(neurons_order, :);
orderedData = orderedData(:, time_order);

% prepare trees for recursion
for treeLevel = 1:length(Trees_all{runningOrder==1})
    row_orderedtree{treeLevel} = Trees_all{runningOrder==1}{treeLevel};
    row_orderedtree{treeLevel}.clustering = Trees_all{runningOrder==1}{treeLevel}.clustering( neurons_order);
end
for treeLevel = 1:length(Trees_all{runningOrder==2})
    col_orderedtree{treeLevel} = Trees_all{runningOrder==2}{treeLevel};
    col_orderedtree{treeLevel}.clustering = Trees_all{runningOrder==2}{treeLevel}.clustering( time_order);
end

figure;
subplot(2,2,1);
imagesc(meanData);title('Orig Data');colormap gray;
xlabel('Time');
set(gca, 'YTick', 1:length(NeuronsLabels));
set(gca, 'YTickLabel',NeuronsLabels);
subplot(2,2,4);
imagesc(orderedData);title('Ordered Data');colormap gray;
xlabel('Time');
set(gca, 'YTick', 1:length(NeuronsLabels));
set(gca, 'YTickLabel',NeuronsLabels(neurons_order));

subplot(2,2,3);
plotTreeWithColors(row_orderedtree, 1:size(meanData,1));
title('Neurons Tree'); colorbar('off')
view(-90,90)
subplot(2,2,2);
plotTreeWithColors(col_orderedtree, 1:size(meanData,2))
title('Time Tree');colorbar('off')

mkNewFolder('D30tilingSols');
files_prefix = 'D30tilingSols\solsD30';

ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
solutionTiling = [];
figure;
clc;
vol_v = setdiff(1e4:-1:64, size(meanData)) ;
getAll2DtilingSols(orderedData, row_orderedtree, col_orderedtree, vol_v, files_prefix, params.verbose);

p.over_rows = true;
p.normalization_type = 'by_std';
orderedDataN = NormalizeData(orderedData, p);
tauvec = [0.1 .5 1 1.5 1.9];
p1.verbose = 2;p1.beta_row = 0;p1.beta_col = 0;
for b = 1:length(tauvec)
for k = 1:length(dir([files_prefix '*']))
  load([files_prefix '_' num2str(k)]);
  p1.tau = tauvec(b);
currErr(k,b) = evalTilingErrTauMeas(orderedDataN, sol, p1);
end
end
% for b = 1:length(betavec)
% for k = 1:length(dir([files_prefix '*']))
%   load([files_prefix '_' num2str(k)]);
%   p1.beta_col = betavec(b);
%   p1.beta_row = betavec(b);
% %   currErr(k,b) = evalTilingErr(orderedDataN, sol, p1);
%   
%   
% %   currErr(k,b) = evalTilingErr_rankCov(orderedDataN, sol, betavec(b), betavec(b));
% end
% [~, minind(b)]=min(currErr(:,b));
% load([files_prefix '_' num2str(minind(b))]);
% figure;imagesc(sol.isbusy);title([num2str(betavec(b)) ' ' num2str(currErr(k,b)) ' ' num2str(max(max(sol.isbusy)))]);
% end
%% plot data and tiling of the other experiments

meanData1 = mean(dataN_all1, 3);
meanData2 = mean(dataN_all2, 3);
meanData3 = mean(dataN_all3, 3);

orderedData1 = meanData1(neurons_order, :);
orderedData1 = orderedData1(:, time_order);
orderedData2 = meanData2(neurons_order, :);
orderedData2 = orderedData2(:, time_order);
orderedData3 = meanData3(neurons_order, :);
orderedData3 = orderedData3(:, time_order);
orderedDataN1 = NormalizeData(orderedData1, p);
orderedDataN2 = NormalizeData(orderedData2, p);
orderedDataN3 = NormalizeData(orderedData3, p);
for b = 1:length(tauvec)
[~, minind(b)]=min(currErr(:,b));
load([files_prefix '_' num2str(minind(b))]);
p1.beta_col = 0;
  p1.beta_row = 0;
  p1.tau = tauvec(b);
[~, meanTiled] = evalTilingErr(orderedData, sol, p1);
% [~, meanTiled1] = evalTilingErr(orderedData1, sol, betavec(b), betavec(b));
% [~, meanTiled2] = evalTilingErr(orderedData2, sol, betavec(b), betavec(b));
% [~, meanTiled3] = evalTilingErr(orderedData3, sol, betavec(b), betavec(b));

figure;plotTiledData(orderedData, meanTiled, NeuronsLabels(neurons_order), sol.isbusy, ['8/12/14 \tau = ' num2str(p1.tau)])
% figure;plotTiledData(orderedData1, meanTiled1, NeuronsLabels(neurons_order), sol.isbusy, '8/15/13')
% figure;plotTiledData(orderedData2, meanTiled2, NeuronsLabels(neurons_order), sol.isbusy, '8/17/14 part 1')
% figure;plotTiledData(orderedData3, meanTiled3, NeuronsLabels(neurons_order), sol.isbusy, '8/17/14 part 2')
% seqMap = (meanTiled -mean(meanTiled(:))) > std(meanTiled(:))/3;
% figure;imagesc(seqMap.*meanTiled);colormap gray;title('D30 - Activation Map');
% set(gca, 'Ytick', 1:length(NeuronsLabels));
% set(gca, 'YtickLabel', NeuronsLabels(neurons_order));

end

% for b = 1:length(betavec)
% [~, minind(b)]=min(currErr(:,b));
% load([files_prefix '_' num2str(minind(b))]);
% p1.beta_col = betavec(b);
%   p1.beta_row = betavec(b);
% [~, meanTiled] = evalTilingErr(orderedData, sol, p1);
% % [~, meanTiled1] = evalTilingErr(orderedData1, sol, betavec(b), betavec(b));
% % [~, meanTiled2] = evalTilingErr(orderedData2, sol, betavec(b), betavec(b));
% % [~, meanTiled3] = evalTilingErr(orderedData3, sol, betavec(b), betavec(b));
% 
% figure;plotTiledData(orderedData, meanTiled, NeuronsLabels(neurons_order), sol.isbusy, '8/12/14')
% % figure;plotTiledData(orderedData1, meanTiled1, NeuronsLabels(neurons_order), sol.isbusy, '8/15/13')
% % figure;plotTiledData(orderedData2, meanTiled2, NeuronsLabels(neurons_order), sol.isbusy, '8/17/14 part 1')
% % figure;plotTiledData(orderedData3, meanTiled3, NeuronsLabels(neurons_order), sol.isbusy, '8/17/14 part 2')
% seqMap = (meanTiled -mean(meanTiled(:))) > std(meanTiled(:))/3;
% figure;imagesc(seqMap.*meanTiled);colormap gray;title('D30 - Activation Map');
% set(gca, 'Ytick', 1:length(NeuronsLabels));
% set(gca, 'YtickLabel', NeuronsLabels(neurons_order));
% 
% end
plotNeuronsClustersWithMartkers(allMat_all);title('8/12/14');
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels(neurons_order));

plotNeuronsClustersWithMartkers(allMat_all1);title('8/15/13');
plotNeuronsClustersWithMartkers(allMat_all2);title('8/17/14 part 1');
plotNeuronsClustersWithMartkers(allMat_all3);title('8/17/14 part 2');

normalized_centroids_all = plotByClustering(allMat_all(:),  ['8/12/14 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all1 = plotByClustering(allMat_all1(:),  ['8/15/13 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all2 = plotByClustering(allMat_all2(:),  ['8/17/14 part 1 Organized By Tree Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)]);%selectedTimeFrams
normalized_centroids_all3 = plotByClustering(allMat_all3(:),  ['8/17/14 part 2 Organized By Tree Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)]);%selectedTimeFrams

% meanvals = flipud(unique(meanTiled(:)));
% for n=1:10
%     [a1, a2] = find(meanTiled == meanvals(n));
%     figure;imagesc(orderedData(min(a1):max(a1), min(a2):max(a2)));
% end
