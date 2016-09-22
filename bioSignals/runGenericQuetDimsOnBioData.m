close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));
experiment = 'D30';
dims4quest = 3;
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
        nt = 360;
        ind2data =3;
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

switch experiment
    case 'D8'
        params.tree{1}.runOnEmbdding = false;
        params.tree{2}.runOnEmbdding = true;
        params.tree{3}.runOnEmbdding = false;
%         selectedTimeFrams = 35:60;
        selectedTimeFrams=1:size(X, 2);
        data = X(:, selectedTimeFrams, :);
        [nr, nt, nT] = size(data);
        toneVec=toneVec(selectedTimeFrams);
        
        neuron_tree_level = 2;
        %% Run Qu.
        params = SetGenericDimsQuestParams(dims4quest, true);
        for ind = 1:dims4quest
            params.emd{ind}.beta = 1;
            
            
        end
        
        params.tree{3}.splitsNum = 2;
            params.tree{3}.treeDepth = inf;
            
          
             params.tree{1}.splitsNum = 9;
            params.tree{1}.treeDepth = inf;
             params.tree{2}.splitsNum = 9;
            params.tree{2}.treeDepth = inf;
            
        if dims4quest == 3
            runningOrder = [3 1 2];
        else
            runningOrder = [1  2];
        end
        params.data.over_rows = false;
        params.data.to_normalize = false;
        params.data.normalization_type = 'by_std';
%         params.init_aff{3}.metric = 'euc';
        % params.init_aff{3}.metric = 'cosine_similarityOnTrials';
        params.init_aff{3}.metric = 'cosine_similarity';
        params.n_iters = 1;
    case 'M2'
        params.tree{1}.runOnEmbdding = true;
        params.tree{2}.runOnEmbdding = false;
        params.tree{3}.runOnEmbdding = false;

        data = X(:, ind2data:3:end, :);
%         selectedTimeFrams = 35:60;
selectedTimeFrams=1:size(data, 2);
        data = data(:, selectedTimeFrams, :);
        [nr, nt, nT] = size(data);
toneVec=toneVec(selectedTimeFrams);
        neuron_tree_level = 2;
        % Run Qu.
        params = SetGenericDimsQuestParams(dims4quest, true);
        for ind = 1:dims4quest
            params.emd{ind}.beta = 0.5;
            params.tree{ind}.splitsNum = 2;
            params.tree{ind}.treeDepth = inf;
            
        end
        if dims4quest == 3
            runningOrder = [1 2 3];
        else
            runningOrder = [1  2];
        end
        params.data.over_rows = false;
        params.data.to_normalize = true;
        params.data.normalization_type = 'by_std';
        % params.init_aff{3}.metric = 'euc';
        % params.init_aff{3}.metric = 'cosine_similarityOnTrials';
        params.init_aff{3}.metric = 'cosine_similarity';
        params.n_iters = 1;
    case 'D30'
        params.tree{1}.runOnEmbdding = true;
        params.tree{2}.runOnEmbdding = false;
        params.tree{3}.runOnEmbdding = false;
        selectedTimeFrams = 115:140;
selectedTimeFrams=1:size(X, 2);
        data =  X(:, selectedTimeFrams, :);
        neuron_tree_level =2;
        %% Run Qu.
        params = SetGenericDimsQuestParams(dims4quest, true);
        for ind = 1:dims4quest
            params.emd{ind}.beta = 1;
            params.tree{ind}.splitsNum = 4;
            params.tree{ind}.treeDepth = 5;
            
        end
        
%         params.tree{3}.splitsNum = 5;
%             params.tree{3}.treeDepth = inf;
            
%             params.tree{1}.splitsNum = 7;
%             params.tree{1}.treeDepth = 3;
%             
%          params.tree{2}.splitsNum = 5;
%             params.tree{2}.treeDepth = 3;
            
        if dims4quest == 3
            runningOrder = [  1 2 3];
        else
            runningOrder = [1  2];
        end
        params.data.over_rows = true;
        params.data.to_normalize = false;
        params.data.normalization_type = 'by_std';
        params.init_aff{3}.metric = 'euc';
        % params.init_aff{3}.metric = 'cosine_similarityOnTrials';
        params.init_aff{3}.metric = 'cosine_similarity';
        params.n_iters = 1;
        
end

params.verbose = 2;

[ Trees, dual_aff, init_aff ] = RunGenericDimsQuestionnaire( params, permute(data,(runningOrder) ) );
% measure the accuracy for neurons level
if strcmp(experiment, 'M2')
    Acc1 = sum((Trees{runningOrder==1}{end-1}.clustering(1:46)==1))/46;
    Acc2 = sum((Trees{runningOrder==1}{end-1}.clustering(47:end)==2))/35;
    figure;plotTreeWithColors(Trees{runningOrder==1}, [ones(1,46) 100*ones(1, 35)]);
    title(['Neurons Tree, Colored By Layer']);
end

figure;
[vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==2}, 0.0), 3);
subplot(2,1,1);
plotEmbeddingWithColors(vecs * vals, 1:size(data, 2), 'Time - ');
subplot(2,1,2);
plotEmbeddingWithColors(vecs * vals, toneVec, 'Time Colored By Tone - ');
if strcmp(experiment, 'M2')
    figure;subplot(2,1,1);
    [vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==1}, 0.0), 3);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Neurons - ');
    subplot(2,1,2);
    plotEmbeddingWithColors(vecs * vals, [ones(1,46) 100*ones(1, 35)], 'Nuerons Colored by Layer - ');
    
else
    figure;
    [vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==1}, 0.0), 3);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 1), 'Nuerons - ');
end
if dims4quest==3
    figure;
    [vecs, vals] = CalcEigs(threshold(dual_aff{runningOrder==3}, 0.0), 2);
    plotEmbeddingWithColors(vecs * vals, 1:size(data, 3), 'Trials - ');
end

%% ordering neurons clusters
clusteringN = Trees{runningOrder==1}{neuron_tree_level}.clustering;
[~, Norder] = sort(clusteringN);
clusteringt = Trees{runningOrder==2}{2}.clustering;
clusteringT = Trees{runningOrder==3}{2}.clustering;

[~, torder] = sort(clusteringt);
[~, Torder] = sort(clusteringT);
clstrsN = unique(clusteringN);
clstrst = unique(clusteringt);

clusteredData = data(Norder, :, :);

figure;imagesc(mean(clusteredData(:, :,:), 3));
y = -.50;x = -.50;
for r = 1:length(clstrsN)
    line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
    
    line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
    y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
end
title([experiment ' - Ordering by Neurons '  'Averaged Over Trials ']);


switch experiment
    %     case 'M2'
    
    case 'D8'
        neuron_tree_level=4;
         getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, :, 1:30), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - Before CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, :, 31:44), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - During CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, :, 45:end), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - After CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       
        selectedTimeFrams = 30:80;
        getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, selectedTimeFrams, 1:30), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - Before CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, selectedTimeFrams, 31:44), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - During CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, selectedTimeFrams, 45:end), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' - After CNO - Organized By Tree Level ' num2str(neuron_tree_level)])
       
        
%         
%         figure;imagesc(mean(clusteredData(:, selectedTimeFrams,1:30), 3));
%         y = -.50;x = -.50;
%         for r = 1:length(clstrsN)
%             line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             
%             line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
%         end
%         title([experiment ' - Ordering by Neurons '  'Averaged Over Trials Before CNO']);
%         figure;imagesc(mean(clusteredData(:, selectedTimeFrams,31:44), 3));
%         y = -.50;x = -.50;
%         for r = 1:length(clstrsN)
%             line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             
%             line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
%         end
%         title([experiment ' - Ordering by Neurons '  'Averaged Over Trials During CNO']);
%         
%         figure;imagesc(mean(clusteredData(:, selectedTimeFrams,45:end), 3));
%         y = -.50;x = -.50;
%         for r = 1:length(clstrsN)
%             line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             
%             line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
%             y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
%         end
%         title([experiment ' - Ordering by Neurons '  'Averaged Over Trials After CNO']);
%         
%         
%         
%         
%         
%         
%         
    case 'D30'
%         selectedTimeFrams = 115:140;
        [X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, 360);
        [X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, 360);
        [X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, 360);
        
        neuron_tree_level = 3;
        getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, :, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' -  8/12/14 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)])
       
        
        getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, data(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels, params.init_aff{runningOrder(1)} , [experiment ' -  8/12/14 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, X1(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels1, params.init_aff{runningOrder(1)} , [experiment ' - 8/15/13 Organized By Tree Of 8/12/14 Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, X2(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels2, params.init_aff{runningOrder(1)} , [experiment ' - 8/17/14 part 1 Organized By Tree Of 8/12/14 Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)])
       getClusteringByTreeAndPlot(Trees{runningOrder==1}{neuron_tree_level}, X3(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels3, params.init_aff{runningOrder(1)} , [experiment ' - 8/17/14 part 2 Organized By Tree Of 8/12/14 Of 8/12/14; Organized By Tree Level ' num2str(neuron_tree_level)])
       
        clusteredData1 = X1(Norder, :, :);
        clusteredData2 = X2(Norder, :, :);
        clusteredData3 = X3(Norder, :, :);
        
        figure;imagesc(mean(clusteredData(:, 115:140,:), 3));
        y = -.50;x = -.50;
        for r = 1:length(clstrsN)
            line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            
            line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
        end
        title(['Ordering by Neurons '  'Averaged Over Trials ']);
        
        figure;imagesc(mean(clusteredData1(:, 115:140,:), 3));
        y = -.50;x = -.50;
        for r = 1:length(clstrsN)
            line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            
            line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
        end
        title(['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
        
        figure;imagesc(mean(clusteredData2(:, 115:140,:), 3));
        y = -.50;x = -.50;
        for r = 1:length(clstrsN)
            line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            
            line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
        end
        title(['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
        
        
        figure;imagesc(mean(clusteredData3(:, 115:140,:), 3));
        y = -.50;x = -.50;
        for r = 1:length(clstrsN)
            line([1 35], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            
            line([1 size(data, 2)], [1 1]*y, 'Color', 'k', 'LineWidth',3);
            y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
        end
        title(['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
        
        
        %
        %     figure;imagesc(mean(clusteredData1, 3));
        %     y = -.50;x = -.50;
        %     for r = 1:length(clstrsN)
        %         line([1 size(data, 2)], [1 1]*y, 'Color', 'k');
        %         y = y + Trees{runningOrder==1}{neuron_tree_level}.folder_sizes(r);
        %     end
        %     for t = 1:length(clstrst)
        %         line([1 1]*x, [1 size(data, 1)],  'Color', 'k');
        %         x = x + Trees{runningOrder==2}{2}.folder_sizes(t);
        %     end
        %     title(['Ordering by Neurons And Time '  'Averaged Over Trials ']);
        %
        
        
        
        % clustering nuerons for every trial
        [~, orderN] = sort(Trees{runningOrder(1)}{neuron_tree_level}.clustering);
        [~, ordert] = sort(Trees{runningOrder(2)}{2}.clustering);
        if dims4quest==3
            
            [~, orderT] = sort(Trees{runningOrder(3)}{2}.clustering);
            clstrsT = unique(Trees{runningOrder(3)}{2}.clustering);
            
        else
            clstrsT = 1:size(data, 3);
        end
        clstrsN = unique(Trees{runningOrder(1)}{neuron_tree_level}.clustering);
        
        clusteredData = data(orderN, :, :);
        % clusteredData = clusteredData(:, ordert, :);
        clusteredData = clusteredData(:, :, orderT);
        
        %     figure;imagesc(mean(clusteredData, 3));
        y = 0;
        allMat = cell(length(clstrsN),1);
        meanMat = zeros(length(clstrsN), size(data, 2));
        for r = 1:length(clstrsN)
            meanMat(r, :) = mean(mean(data(Trees{runningOrder(1)}{neuron_tree_level}.clustering==clstrsN(r), :, :), 1), 3);
            allMat{r} = mean(data(Trees{runningOrder(1)}{neuron_tree_level}.clustering==clstrsN(r), :, :), 3);
        end
        
        % plotByClustering(meanMat, allMat, '', 30:60)
        % plotByClustering(meanMat,  '', 110:140);
        % plotByClustering(allMat, '', 110:140);
        plotByClustering(allMat, '', 1:31);
end
%     case 'D30'
        [X1, ~, NeuronsLabels1] = loadNeuronsData(datapth, {'8_15_13_1-35'}, 360);
        [X2, ~, NeuronsLabels2] = loadNeuronsData(datapth, {'8_17_14_1-45'}, 360);
        [X3, ~, NeuronsLabels3] = loadNeuronsData(datapth, {'8_17_14_46-80'}, 360);

        [meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(Trees{runningOrder(1)}{neuron_tree_level}, data, NeuronsLabels, NeuronsLabels);
        [meanMat1, allMat1, meanMatAlltrials1] = getCentroidsByTree(Trees{runningOrder(1)}{neuron_tree_level}, X1(:, :, :), NeuronsLabels, NeuronsLabels1);
        [meanMat2, allMat2, meanMatAlltrials2] = getCentroidsByTree(Trees{runningOrder(1)}{neuron_tree_level}, X2(:, :, :), NeuronsLabels, NeuronsLabels2);
        [meanMat3, allMat3, meanMatAlltrials3] = getCentroidsByTree(Trees{runningOrder(1)}{neuron_tree_level}, X3(:, :, :), NeuronsLabels, NeuronsLabels3);
        [aff_mat1  ] = CalcInitAff( meanMat.', params.init_aff{runningOrder(1)} );
        [vecs, vals] = CalcEigs(threshold(aff_mat1, 0.0), 2);%
        [ row_order1 ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );
        % figure;plotEmbeddingWithColors(vecs(row_order,:) * vals, 1:11, 'Nuerons Embedding');
%         plotByClustering(meanMat(row_order1, :),  ['8/12/14 Organized By Tree Level ' num2str(neuron_tree_level)]);
%         plotByClustering(meanMat1(row_order1, :),  ['8/15/13 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
%         plotByClustering(meanMat2(row_order1, :),  ['8/17/14 part 1 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
%         plotByClustering(meanMat3(row_order1, :), ['8/17/14 part 2 Organized By Tree Of 8/12/14; Level ' num2str(neuron_tree_level)]);
%
        plotByClustering(allMat(row_order1),  ['8/12/14 Organized By Tree Level ' num2str(neuron_tree_level)]);
        plotByClustering(allMat1(row_order1),  ['8/12/14 Organized By Tree Level ' num2str(neuron_tree_level)]);
        plotByClustering(allMat2(row_order1),  ['8/12/14 Organized By Tree Level ' num2str(neuron_tree_level)]);
        plotByClustering(allMat3(row_order1),  ['8/12/14 Organized By Tree Level ' num2str(neuron_tree_level)]);
%
%
%         end
