close all;
clear all;
clc;
dbstop if error;

files = {  '8_6_14_1-20_control' '8_4_14_1-25_control'};%'8_6_14_21-60_cno' 
rng(73631);
experiment = 'D8';
dims4quest = 3;

datapth = '..\..\..\datasets\biomed\D8';
files = { '8_6_14_1-20_control' '8_6_14_21-60_cno'};
nt = 120;
toneVec = [ones(1,40) 100*ones(1, 80)];
[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);

selectedTimeFrams=2:size(X, 2);
data = X(:, selectedTimeFrams, :);
toneVec=toneVec(selectedTimeFrams);
istopdown = true;

params = SetGenericDimsQuestParams(dims4quest, istopdown);
    runningOrder = [3 1 2];
    params.data.over_rows = false;

    for ind = 1:dims4quest
        params.emd{ind}.beta = .1;
        params.tree{ind}.splitsNum = 2;
        params.tree{ind}.treeDepth = 6;
        
    end
    params.data.over_rows = true;
    params.n_iters = 1;

    params.init_aff{3}.thresh = .1;
    params.tree{1}.runOnEmbdding = true;
params.tree{2}.runOnEmbdding = true;
params.tree{3}.runOnEmbdding = true;
params.data.to_normalize = true;
params.data.normalization_type = 'by_std';
params.verbose = 2;
[ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params, permute(data,(runningOrder) ) );

figure;
subplot(2,1,1);
plotTreeWithColors(Trees{runningOrder == 2}, toneVec);
title('Time Tree; Colored By Tone');

[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==2}, 0.1), 4);
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Time Colored By Tone - ');view(22, 32);

figure;
subplot(2,1,1);
plotTreeWithColors(Trees{runningOrder == 3}, [ones(25,1); 2*ones(0, 1); 3*ones(20,1); 4*ones(15,1)]);
title('Trials Tree; Colored By CNO State (Before, During & After)');
[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==3}, 0.0), 3);
subplot(2,1,2);plotEmbeddingWithColors(vecs * vals, [ones(25,1); 2*ones(0, 1); 3*ones(20,1); 4*ones(15,1)], 'Trials - ');

%% Order by tree
dataN = NormalizeData(permute(data, runningOrder), params.data);
[~, ic] = sort(runningOrder);
dataN = permute(dataN, ic);

neuron_tree_level = 2;
[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(Trees{runningOrder==1}{neuron_tree_level}, dataN(:,:,1:10), NeuronsLabels, NeuronsLabels);
affine = feval(params.init_aff{runningOrder==1}.initAffineFun, permute(meanMatAlltrials, runningOrder), params.init_aff{1});
[vecs, vals] = CalcEigs(affine, 5);
classes_order = OrganizeDiffusion3DbyOneDim( meanMatAlltrials, vecs*vals );

normalized_centroids = plotByClustering(allMat(classes_order),  '');%selectedTimeFrams
      

[~, neurons_order] = sort(Trees{runningOrder==1}{2}.clustering);
[~, time_order] = sort(Trees{runningOrder==2}{2}.clustering);
orderedData = dataN(neurons_order, :, 1);
orderedData = orderedData(:, time_order, 1);

% prepare trees for recursion
for treeLevel = 1:length(Trees{runningOrder==1})
    row_orderedtree{treeLevel} = Trees{runningOrder==1}{treeLevel};
    row_orderedtree{treeLevel}.clustering = Trees{runningOrder==1}{treeLevel}.clustering( neurons_order);
end
for treeLevel = 1:length(Trees{runningOrder==2})
    col_orderedtree{treeLevel} = Trees{runningOrder==2}{treeLevel};
    col_orderedtree{treeLevel}.clustering = Trees{runningOrder==2}{treeLevel}.clustering( time_order);
end

mkNewFolder('D8tilingSols');
files_prefix = 'D8tilingSols\solsD30_';

ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
solutionTiling = [];
figure;
clc;
l = 1;
vol_v = {setdiff([  1000:-1:100  ], size(orderedData))};
for vol_i = 1:length(vol_v)
    [minCurrErr, currSolutionTiling] = loopTiling2D(orderedData, row_orderedtree, col_orderedtree, vol_v{vol_i}, files_prefix);

%     [minCurrErr, tilingCurrRes, currSolutionTiling] = recursiveTiling2D(orderedData(:,:,1), row_orderedtree, col_orderedtree, vol_v(vol_i), ind2data, tiling, solutionTiling, Inf);
    if minCurrErr < inf
        volumeRes(l) = vol_v(vol_i);
        minErr(l) = minCurrErr;
        solutionTilingRes(l) = currSolutionTiling;
        l = l + 1;
    end
end

files = 'solsD8_';
betavec = -1:.2:1;
for b = 1:length(betavec)
%     figure;
    
for k = 1:15
%     subplot(5,3,k)
  load([files num2str(k)]);
  currErr(k,b) = evalTilingErr(orderedData, sol, betavec(b), betavec(b));
%   imagesc(sol.isbusy);title([num2str(betavec(b)) ' ' num2str(currErr(k,b)) ' ' num2str(max(max(sol.isbusy)))]);
  
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




