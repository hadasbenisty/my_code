function neuronsClusteringConsistencyAllAnimals(animalName, figspath, th )


dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));
addpath(genpath('../../grangerCausality/'));
datapth = ['../../../datasets/biomed/' animalName ];

nt = 360;
switch animalName
    case 'D13'
        files = {'7_21_14_1-30' '7_1_14_1-16',  '6_28_14_1-21' };
        trialsVec = [ones(20,1); ones(16,1)*2; ones(30,1)*3];
        toneTimeFrames = 110:150;%115:140;

    case 'D30'
        files = {'8_12_14_1-40' '8_15_13_1-35' '8_17_14_1-45' '8_17_14_46-80'};
        trialsVec = [ones(40,1); ones(35,1)*2; ones(45,1)*3; ones(35,1)*4];
        toneTimeFrames = 110:150;%115:140;

    case 'M2'
        files = {'4_4_14'};
    case 'D8'
        files = { '8_6_14_1-20_control' '8_6_14_21-60_cno' };
        tileAll = false;
        nt = 120;
        toneTimeFrames = 30:60;
end

[expDates, namesNums] = filename2date(files);


toneVec = [ones(1,120) 100*ones(1, 240)];
%% Load Data
X = cell(length(files), 1); NeuronsLabels = cell(length(files), 1);
for n = 1:length(files)
    [X{n}, ~, NeuronsLabels{n}] = loadNeuronsData(datapth, files(n), nt);
end

%% Set Params
% selectedTimeFrams = 115:140;

dims4quest=3;
params = SetGenericDimsQuestParams(dims4quest, true);
for ind = 1:dims4quest
    %             params.init_aff{ind}.metric = 'cosine_similarity';
    params.init_aff{ind}.metric = 'euc';
    params.tree{ind}.runOnEmbdding = true;
    params.tree{ind}.splitsNum=[2];
    params.tree{ind}.treeDepth = 8;%4;
end
params.tree{2}.eigs_num=5;
params.tree{2}.splitsNum = [10 2 2 2 2 2 2 2];%9
runningOrder = [  3 2 1];


params.data.over_rows = ~true;
params.data.to_normalize = true;%false
params.data.normalization_type = 'by_std';

X_N = cell(length(files), 1);
for n = 1:length(X)
    X_N{n} = normalizeDataWithOrdering(X{n}, runningOrder, params.data);
end
params.data.to_normalize = false;

params.n_iters = 2;
params.verbose = 2;


%% Run Qu.

Trees = cell(length(files), 1); dual_aff = cell(length(files), 1); init_aff = cell(length(files), 1);
for n = 1
    [ Trees{n}, dual_aff{n}, init_aff{n} ] = RunGenericDimsQuestionnaire( params, permute(X_N{n},(runningOrder) ) );
end


%% Get Sequence By Tree Ordering

neuron_tree_level = 3;
meanData = cell(length(X_N), 1); meanMat = cell(length(X_N), 1); allMat = cell(length(X_N), 1); meanMatAlltrials = cell(length(X_N), 1);
for n = 1:length(X_N)
    meanData{n} = mean(X_N{n}, 3);
    [meanMat{n}, allMat{n}, meanMatAlltrials{n}] = getCentroidsByTree(Trees{1}{runningOrder==1}{neuron_tree_level}, meanData{n}, NeuronsLabels{1}, NeuronsLabels{n});
end

affine{1} = feval(params.init_aff{runningOrder==1}.initAffineFun, meanMatAlltrials{1}.', params.init_aff{1}, find(runningOrder==2));
[vecs, vals] = CalcEigs(affine{1}, 9);
classes_order{1} = OrganizeDiffusion3DbyOneDim( meanMatAlltrials{1}, vecs*vals );
normalized_centroids = cell(length(X_N), 1);
n=1;
mkNewFolder(fullfile(figspath, animalName));

normalized_centroids{n} = plotByClustering(allMat{n}(classes_order{1}),  [animalName ' - ' expDates{n} ]);%selectedTimeFrams
mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'centroids']));
xlim([toneTimeFrames(1), toneTimeFrames(end)]);
mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'centroidsZoom']));
ca;
plotOnsetByClustering(meanMat{n}(classes_order{1}, :),  [animalName ' - ' expDates{n} ], 1:size(normalized_centroids{n}, 2), th,true, 3);
mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'activation']));
f=gcf;a = get(f,'Children');
xlim(a(1),[toneTimeFrames(1), toneTimeFrames(end)]);xlim(a(2),[toneTimeFrames(1), toneTimeFrames(end)]);
mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'activationZoom']));
ca;
% plotNeuronsClustersWithMartkers(allMat{n});title([animalName ' - ' expDates{n} ]);

for n = 2:length(X_N)
    normalized_centroids{n} = plotByClustering(allMat{n}(classes_order{1}),  [animalName ' - ' expDates{n} ' Organized By Tree Of ' files{1}]);%selectedTimeFrams
    mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'centroids']));
    xlim([toneTimeFrames(1), toneTimeFrames(end)]);
    mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'centroidsZoom']));
    ca;
    plotOnsetByClustering(meanMat{n}(classes_order{1}, :),  [animalName ' - ' expDates{n} ' Organized By Tree Of ' files{1}], 1:size(normalized_centroids{n}, 2), th,true, 3);
    mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'activation']));
    f=gcf;a = get(f,'Children');
    xlim(a(1),[toneTimeFrames(1), toneTimeFrames(end)]);xlim(a(2),[toneTimeFrames(1), toneTimeFrames(end)]);
    mysave(gcf, fullfile(figspath, animalName, [namesNums{n} 'activationZoom']));
    ca;
     
   
    %     plotNeuronsClustersWithMartkers(allMat{n});title([animalName ' - ' expDates{n} ' Organized By Tree Of ' files{1}]);
end






