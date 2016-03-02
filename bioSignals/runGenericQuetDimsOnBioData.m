close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../gen_utils'));
addpath(genpath('../tiling'));
experiment = 'M2';
switch experiment
    case 'D30'
        datapth = '..\..\..\datasets\biomed\D30';
        nt = 360;
        toneVec = [ones(1,120) 100*ones(1, 240)];
        files = { '8_12_14_1-40' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'
    case 'D8'
        datapth = '..\..\..\datasets\biomed\D8';
        files = { '8_6_14_1-20_control' '8_6_14_21-60_cno'};
        nt = 120;
        toneVec = [ones(1,40) 100*ones(1, 80)];

    case 'M2'
        datapth = '..\..\..\datasets\biomed\M2';
        nt = 120;
        toneVec = [ones(1,40) 100*ones(1, 80)];
        files = { '4_4_14' };% '8_15_13_1-35' '8_12_14_1-40' '8_17_14_1-45'  '8_17_14_46-80'  '8_15_13_1-35' '8_12_14_1-40'

end
rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
%% Load Data


[X, expLabel, NeuronsLabels] = loadNeuronsData(datapth, files, nt);

% selectedTimeFrams = 110:140;
selectedTimeFrams=1:size(X, 2);
Xsel = X(:, selectedTimeFrams, :);
[nr, nt, nT] = size(Xsel);
data = Xsel;
%% Run Qu. 2D
params = SetGenericDimsQuestParams(2, true);
for ind = 1:2
    params.emd{ind}.beta = 0.5;
    params.tree{ind}.treeDepth = 8;
    
end
for ciT = 1:size(data, 3);
    params.verbose = 0;
    [ Trees2D, dual_aff2D, init_aff2D ] = RunGenericDimsQuestionnaire( params, permute(data(:,:,ciT), [2 1]) );

    
    clusteringN = Trees2D{2}{2}.clustering;
    [~, Norder] = sort(clusteringN);
    clusteringt = Trees2D{1}{2}.clustering;
    [~, torder] = sort(clusteringt);
    
    clstrsN = unique(clusteringN);
    clstrst = unique(clusteringt);
    
    clusteredData2D(:,:,ciT) = data(Norder, :, ciT);
    clusteredData2D1(:,:,ciT) = clusteredData2D(:, torder, ciT);



    figure;imagesc(clusteredData2D(:,:,ciT));
    y = -.50;x = -.50;
    for r = 1:length(clstrsN)
        line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
        y = y + Trees2D{2}{2}.folder_sizes(r);
    end
    title(['Trial no. ' num2str(ciT) ' Ordering by Neurons']);
    
    
    
%     figure;imagesc(clusteredData2D1(:,:,ciT));
%     y = -.50;x = -.50;
%     for r = 1:length(clstrsN)
%         line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
%         y = y + Trees2D{2}{2}.folder_sizes(r);
%     end
%     for t = 1:length(clstrst)
%         line([1 1]*x, [1 size(data, 1)],  'Color', 'k');
%         x = x + Trees2D{1}{2}.folder_sizes(t);
%     end
%     title(['Trial no. ' num2str(ciT) ' Ordering by Neurons And Time']);
    
    
    
end
%% Run Qu. 3D
params = SetGenericDimsQuestParams(3, true);
for ind = 1:3
    params.emd{ind}.beta = 0.5;
    params.tree{ind}.treeDepth = 8;
    
end
[ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params, permute(data, [2 1 3]) );
for ciT = 1:size(data, 3);

clusteringN = Trees{2}{2}.clustering;
    [~, Norder] = sort(clusteringN);
    clusteringt = Trees{1}{2}.clustering;
    [~, torder] = sort(clusteringt);
    
    clstrsN = unique(clusteringN);
    clstrst = unique(clusteringt);
    
    clusteredData(:,:,ciT) = data(Norder, :, ciT);
    clusteredData1(:,:,ciT) = clusteredData(:, torder, ciT);



    figure;imagesc(clusteredData(:,:,ciT));
    y = -.50;x = -.50;
    for r = 1:length(clstrsN)
        line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
        y = y + Trees{2}{2}.folder_sizes(r);
    end
    title(['Trial no. ' num2str(ciT) ' Ordering by Neurons']);
    
    
    
%     figure;imagesc(clusteredData2D1(:,:,ciT));
%     y = -.50;x = -.50;
%     for r = 1:length(clstrsN)
%         line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
%         y = y + Trees2D{2}{2}.folder_sizes(r);
%     end
%     for t = 1:length(clstrst)
%         line([1 1]*x, [1 size(data, 1)],  'Color', 'k');
%         x = x + Trees2D{1}{2}.folder_sizes(t);
%     end
%     title(['Trial no. ' num2str(ciT) ' Ordering by Neurons And Time']);
    
end
figure;
[vecs, vals] = CalcEigs(threshold(dual_aff{1}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time Embedding');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Time Embedding');

figure;
[vecs, vals] = CalcEigs(threshold(dual_aff{2}, 0.0), 3);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons Embedding');
figure;
[vecs, vals] = CalcEigs(threshold(dual_aff{3}, 0.0), 4);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials Embedding');


% clustering nuerons for every trial
clusteringN = Trees{1}{2}.clustering;
clusteringT = Trees{3}{2}.clustering;

clstrsN = unique(clusteringN);
clstrsT = unique(clusteringT);

clusteredData = data(clusteringN, :, :);
clusteredData = clusteredData(:, :, clusteringT);
for ciT = 1:size(data, 3);
    figure;imagesc(clusteredData(:,:,ciT));
    y = 0;
    for r = 1:length(clstrsN)
        line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
        y = y + Trees{1}{2}.folder_sizes(r);
    end
    title(['Trial no. ' num2str(ciT)]);
end
for ciT = 1:size(data, 3);
    for ciN = 1:length(clstrsN);
        neurons = find(clusteringN == clstrsN(ciN));
        mat = data(neurons, :, ciT);
        clusteredData(1:length(neurons), :, ciT) = mat;
    end
end
% clustering trials and nuerons
clusteringT = Trees{3}{2}.clustering;
clusteringN = Trees{1}{2}.clustering;
clstrsT = unique(clusteringT);
clstrsN = unique(clusteringN);
currData = cell(length(clstrsT), length(clstrsN));

for ciT = 1:length(clstrsT);
    trials = find(clusteringT == clstrsT(ciT));
    for ciN = 1:length(clstrsN);
        neurons = find(clusteringN == clstrsN(ciN));
        mat = data(neurons, :, trials);
        
        for trialsi = 1:length(trials)
            currData{ciT, ciN} = cat(1, currData{ciT, ciN}, mat(:,:,trialsi));
        end
    end
end
figure;
l = 1;
for n = 1:size(currData,1)
    for m = 1:size(currData,2)
        subplot(size(currData,1), size(currData,2), l);
        imagesc(currData{n, m});
        l=l + 1;
    end
end
% [X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, 360);
% [X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, 360);
% [X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, 360);
% X1sel = X1(:, selectedTimeFrams, :);
% X2sel = X2(:, selectedTimeFrams, :);
% X3sel = X3(:, selectedTimeFrams, :);
%
% tree_level = 2;
% [meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(row_tree{tree_level}, Xsel, NeuronsLabels, NeuronsLabels);
%
% figure;
% clusters = unique(trial_tree{tree_level}.clustering);
% R = ceil(sqrt(length(clusters)));
% for ci = 1:length(clusters)
%     subplot(R,R,ci);
%    inds = find(trial_tree{tree_level}.clustering == clusters(ci));
%    meanMatAlltrialsT(:, :, ci) = mean(meanMatAlltrials(:, :, inds), 3);
%    imagesc(meanMatAlltrialsT(:,:,ci));
% end
% [meanMat1, allMat1, meanMatAlltrials1] = getCentroidsByTree(row_tree{tree_level}, X1sel, NeuronsLabels, NeuronsLabels1);
% [meanMat2, allMat2, meanMatAlltrials2] = getCentroidsByTree(row_tree{tree_level}, X2sel, NeuronsLabels, NeuronsLabels2);
% [meanMat3, allMat3, meanMatAlltrials3] = getCentroidsByTree(row_tree{tree_level}, X3sel, NeuronsLabels, NeuronsLabels3);
% figure;
%
%
% figure;
% for T = 1:size(meanMatAlltrials, 3)
%     subplot(2,2,1);    imagesc(permute(meanMatAlltrials(T,:,:),[3 2 1]))
%     subplot(2,2,2);    imagesc(permute(meanMatAlltrials1(T,:,:),[3 2 1]))
%     subplot(2,2,3);    imagesc(permute(meanMatAlltrials2(T,:,:),[3 2 1]))
%     subplot(2,2,4);    imagesc(permute(meanMatAlltrials3(T,:,:),[3 2 1]))
%     pause;
% end
% [ aff_mat ] = CalcEmdAffOnTreeLevels( meanMat.', col_tree, tree_level, params.row_emd);
% [vecs, vals] = CalcEigs(threshold(aff_mat, 0.0), 3);
% [ row_order ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% % figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
% % showing that the initial metric is  as good as the EMD
% %
% [ aff_mat1,  ] = CalcInitAff( meanMat.', params.init_aff_row );
% [vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);%
% [ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
% % figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
% plotByClustering(meanMat(row_order1, :), allMat(row_order1), ['8/12/14 Organized By Tree Level ' num2str(tree_level)]);
% plotByClustering(meanMat1(row_order1, :), allMat1(row_order1), ['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
% plotByClustering(meanMat2(row_order1, :), allMat2(row_order1), ['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
% plotByClustering(meanMat3(row_order1, :), allMat3(row_order1), ['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(tree_level)]);
%
