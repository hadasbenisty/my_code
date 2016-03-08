clc;
clear all;
close all hidden;
close all;
dbstop if error;
addpath(('../../SvdKmeansCluster/'));
addpath('../../3D_Questionnaire/Questionnaire/');
addpath('../../PCAclustering/');
addpath('../../Matlab_Utils/');


params(1).k = 8;
params(1).embedded = true;
params(1).threshold = 0.4;
params(1).title = 'Neurons';
params(1).verbose = 1;

params(2).k = 5;
params(2).embedded = false;
params(2).threshold = 0.05;
params(2).title = 'Time';
params(2).verbose = 1;

params(3).k = linspace(10,50,5);
params(3).embedded = false;
params(3).threshold = 0.7;
params(3).title = 'Trials';
params(3).verbose = 1;


filePrefix = '../../../datasets\biomed\D30/8_12_14_1-40'; % 0.5
[matrix, ~, NeuronsLabels] = loadNeuronsData('../../../datasets\biomed\D30/', {'8_12_14_1-40'}, 360);
load(strcat(filePrefix,'sensors_dat_titles.mat'));
load('../../svdKmeansCluster/D30/commonIndex')
[sizeRows,sizeCols,maxTime] = size(matrix);
shuffleRows = 1:sizeRows;
shuffleCols = 1:sizeCols;
shuffleTime = 1:maxTime;
data = matrix;
dataOrig = data;
%% search parameters
k_t_vec = [2 3 4 5 6 7];
k_r_vec = [5 6 7 8 9 10];
l = 1;

for threshold = 0.4:0.05:0.65
    
    for k_t = k_t_vec
        params(2).k = k_t;
        for k_r = k_r_vec
            params(1).k = k_r;
            params(1).threshold = threshold;
            [class, vectors, affinity] = svdClustering3D(data, params);
            if length(unique(class{1})) == 9
%                 disp([k_t, k_r]);
%                 figure;
%                 [row_vecs, row_vals] = CalcEigs(affinity{1}, 4);
%                 embedding = row_vecs*row_vals;
%                 %             subplot(length(k_t_vec), length(k_r_vec), l);
%                 plotEmbeddingWithColors(embedding, class{1}, num2str(l));
%                 colorbar off;
                disp(['l = ' num2str(l) ' Neurons - k_t = ' num2str(k_t) ' k_r = ' num2str(k_r) ' threshold = ' num2str(threshold)]);
%                 xlim([-0.1 0.3])
%                 xlim([-0.05 0.1])
%                 view([-115 10]);
%                 drawnow;
                dummyTree.clustering = class{1};
                [meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(dummyTree, data, NeuronsLabels, NeuronsLabels);
%                 plotByClustering(allMat,  num2str(l));drawnow;
                plotByClustering(meanMat,  num2str(l));drawnow;                
                l = l + 1;
            end
        end
    end
end

%% Run on selected parameters
% 1  k_t = 3 k_r = 9 threshold = 0.4 not good
% 3  k_t = 5 k_r = 7 threshold = 0.4 best
% 5  k_t = 7 k_r = 8 threshold = 0.4 not good
% 7  k_t = 3 k_r = 8 threshold = 0.45 not good
% 11 k_t = 7 k_r = 9 threshold = 0.45 not good
% 15 k_t = 7 k_r = 9 threshold = 0.5 maybe
params(2).k = 5;
params(1).k = 7;
params(1).threshold = 0.4;
[class, vectors, affinity] = svdClustering3D(data, params);
dummyTree.clustering = class{1};
[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(dummyTree, data, NeuronsLabels, NeuronsLabels);
% plotByClustering(allMat,  ['8/12/14 Organized By SVD Clustering']);


[X1, ~, NeuronsLabels1] = loadNeuronsData('../../../datasets\biomed\D30/', {'8_15_13_1-35'}, 360);
[X2, ~, NeuronsLabels2] = loadNeuronsData('../../../datasets\biomed\D30/', {'8_17_14_1-45'}, 360);
[X3, ~, NeuronsLabels3] = loadNeuronsData('../../../datasets\biomed\D30/', {'8_17_14_46-80'}, 360);

[meanMat1, allMat1, meanMatAlltrials1] = getCentroidsByTree(dummyTree, X1(:, :, :), NeuronsLabels, NeuronsLabels1);
[meanMat2, allMat2, meanMatAlltrials2] = getCentroidsByTree(dummyTree, X2(:, :, :), NeuronsLabels, NeuronsLabels2);
[meanMat3, allMat3, meanMatAlltrials3] = getCentroidsByTree(dummyTree, X3(:, :, :), NeuronsLabels, NeuronsLabels3);

plotByClustering(meanMat,  ['8/12/14 Organized By SVD Clustering']);
plotByClustering(meanMat1,  ['8/15/13 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(meanMat2,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(meanMat3,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);


plotByClustering(allMat,  ['8/15/13 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat,  ['8/12/14 Organized By SVD Clustering']);


% filePrefix = strcat('D30/ClassResults2',filePrefix(4:end)); %0.5

% save(strcat(filePrefix,'affinity'),'affinityRows','classRows');
% save(strcat(filePrefix,'affinityCols'),'affinityCols','classCols');

% [ err_rate,row_order,col_order,trial_order ] = OrganizeData3D( dataOrig, data, affinityRows, affinityCols{1},...
%  affinityTime, shuffleRows, shuffleCols, shuffleTime, 4, 4, 4 );


