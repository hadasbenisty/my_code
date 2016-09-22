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
selectedTimeFrams = 115:140;
data_all = matrix;
data_tone = matrix(:, selectedTimeFrams, :);

%% search parameters
k_t_vec = 5;%[2 3 4 5 6 7];
k_r_vec = 8;%[5 6 7 8 9 10];
l = 1;
if 0
    for threshold = 0.4:0.05:0.65
        
        for k_t = k_t_vec
            params(2).k = k_t;
            for k_r = k_r_vec
                params(1).k = k_r;
                params(1).threshold = threshold;
                [class, vectors, affinity] = svdClustering3D(data_all, params);
                
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
                [meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(dummyTree, data_all, NeuronsLabels, NeuronsLabels);
                %                 plotByClustering(allMat,  num2str(l));drawnow;
                plotByClustering(meanMat,  num2str(l));drawnow;
                l = l + 1;
                
            end
        end
    end
end
%% Run on selected parameters

params(2).k = 5;
params(1).k = 8;
params(1).threshold = 0.4;
[class_all, vectors_all, affinity_all] = svdClustering3D(data_all, params);

params(2).k = 5;
params(1).k = 8;
params(1).threshold = 0.4;
params(3).k = [10, 20];
[class_tone vectors_tone, affinity_tone] = svdClustering3D(data_tone, params);

dummyTree.clustering = class_all{1};
[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(dummyTree, data_all, NeuronsLabels, NeuronsLabels);


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


plotByClustering(allMat,  ['8/12/14 Organized By SVD Clustering']);
plotByClustering(allMat1,  ['8/15/13 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat2,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat3,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);

% zoom into the tone

plotByClustering(meanMat,  ['8/12/14 Organized By SVD Clustering'], selectedTimeFrams);
plotByClustering(meanMat1,  ['8/15/13 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);
plotByClustering(meanMat2,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);
plotByClustering(meanMat3,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);
th_activation=.7;
plotOnsetByClustering(meanMat,  ['8/12/14 Organized By SVD Clustering'], 1:360, th_activation);xlim([100 220]);
% plotOnsetByClustering(allMat,  ['8/12/14 Organized By SVD Clustering'], 1:360, th_activation);xlim([115 140]);

plotOnsetByClustering(meanMat1,  ['8/15/13 Organized By SVD Clustering Of 8/12/14'],  1:360, th_activation);xlim([100 220]);
plotOnsetByClustering(meanMat2,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14'],  1:360, th_activation);xlim([100 220]);
plotOnsetByClustering(meanMat3,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14'],  1:360, th_activation);xlim([100 220]);

plotByClustering(allMat,  ['8/12/14 Organized By SVD Clustering'], selectedTimeFrams);
plotByClustering(allMat1,  ['8/15/13 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);
plotByClustering(allMat2,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);
plotByClustering(allMat3,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14'], selectedTimeFrams);


dummyTree.clustering = class_tone{1};
[meanMat_tone, allMat_tone, meanMatAlltrials] = getCentroidsByTree(dummyTree, data_tone, NeuronsLabels, NeuronsLabels);
[meanMat1_tone, allMat1_tone, meanMatAlltrials1] = getCentroidsByTree(dummyTree, X1(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels1);
[meanMat2_tone, allMat2_tone, meanMatAlltrials2] = getCentroidsByTree(dummyTree, X2(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels2);
[meanMat3_tone, allMat3_tone, meanMatAlltrials3] = getCentroidsByTree(dummyTree, X3(:, selectedTimeFrams, :), NeuronsLabels, NeuronsLabels3);


plotByClustering(meanMat_tone,  ['8/12/14 Organized By SVD Clustering']);
plotByClustering(meanMat1_tone,  ['8/15/13 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(meanMat2_tone,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(meanMat3_tone,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);

plotByClustering(allMat_tone,  ['8/12/14 Organized By SVD Clustering']);
plotByClustering(allMat1_tone,  ['8/15/13 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat2_tone,  ['8/17/14 Part 1 Organized By SVD Clustering Of 8/12/14']);
plotByClustering(allMat3_tone,  ['8/17/14 Part 2 Organized By SVD Clustering Of 8/12/14']);

