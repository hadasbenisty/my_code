clc;
clear all;
close all hidden;
close all;
dbstop if error;
addpath(genpath('../../3D_Questionnaire\Questionnaire'))
addpath('../../SvdKmeansCluster');
embedded   = false;
coupling = true;
randData = false;
plotFlag = false;
plotFlagKmeans = false;
useFit = false;
run_on_suffled_data = true;
dims4quest = 3;
runningOrder = [3 2 1  ];%3 1 2

load('../../../datasets/Animation\workspaceAnimation2');

if run_on_suffled_data
    X = double(data);
    rows_perm = shuffleRows;
    col_perm = shuffleCols;
    trial_perm = shuffleTime;
else
    X = double(dataOrig);
    rows_perm = 1:length(shuffleRows);
    col_perm = 1:length(shuffleCols);
    trial_perm = 1:length(shuffleTime);
end
[sizeRows,sizeCols,maxTime] = size(X);


addpath(genpath('../../PCAclustering/'));
addpath('../../Matlab_Utils/');

%% Ordering by initial metrics;
% params  = SetGenericQuestParamsAnimation;
% 
% col_init_aff = feval(params.init_aff_col.initAffineFun, X, params.init_aff_col);
% trial_init_aff = feval(params.init_aff_col.initAffineFun, permute(X, [2 1 3]), params.init_aff_col);
% row_init_aff = feval(params.init_aff_col.initAffineFun, permute(X, [1 3 2]), params.init_aff_col);
% [ err_rate2,organized_data2,row_order2,col_order2,orderTime2 ] = OrganizeData3D( double(dataOrig), X, trial_init_aff, col_init_aff,...
% row_init_aff, rows_perm, col_perm, trial_perm, 50, 50, 50 );
% organized_data_initmetrics = X(:,col_order2, :); 
% organized_data_initmetrics = organized_data_initmetrics(row_order2,:, :); 
% organized_data_initmetrics = organized_data_initmetrics(:,:, orderTime2);
% implay(uint8(organized_data_initmetrics));

%% run BU trees 
params = SetGenericDimsQuestParams(3, true);
for ind = 1:dims4quest
    params.emd{ind}.beta = 0.0;
    
    
end
% [ TreesBU, dual_affBU, init_affBU ] = RunGenericDimsQuestionnaire( params, permute(data,(runningOrder) ) );
% [ err_rateBU,organized_dataBU,row_orderBU,col_orderBU,orderTimeBU ] = OrganizeData3D( double(dataOrig), double(X), dual_affBU{1}, dual_affBU{2},...
% dual_affBU{3}, rows_perm, col_perm, trial_perm, 50, 50, 50 );
% 
% organized_data_BUtrees = X(:,col_orderBU, :); 
% organized_data_BUtrees = organized_data_BUtrees(row_orderBU,:, :); 
% organized_data_BUtrees = organized_data_BUtrees(:,:, orderTimeBU); 
% implay(uint8(organized_data_BUtrees));

%% run TD trees with SVD clustering
params = SetGenericDimsQuestParams(3, true);
for n = 1:3
    params.tree{n}.splitsNum = 9;
        params.emd{ind}.beta = 0;
    params.tree{ind}.treeDepth = 4;
    params.tree{ind}.runOnEmbdding = false;
end
params.data.over_rows = false;
params.data.to_normalize = false;
params.data.normalization_type = 'by_std';
[ TreesTD, dual_affTD, init_affTD ] = RunGenericDimsQuestionnaire( params, permute(X,(runningOrder) ) );

[ err_rateTD,organized_dataTD,row_orderTD,col_orderTD,orderTimeTD ] = OrganizeData3D( double(dataOrig), double(X), dual_affTD{runningOrder==1}, dual_affTD{runningOrder==2},...
dual_affTD{runningOrder==3}, rows_perm, col_perm, trial_perm, 12, 12, 12 );

organized_data_TDtrees = X(:,col_orderTD, :); 
organized_data_TDtrees = organized_data_TDtrees(row_orderTD,:, :); 
organized_data_TDtrees = organized_data_TDtrees(:,:, fliplr(orderTimeTD)); 
implay(uint8(organized_data_TDtrees));



% columns - apply svdClass for different scales
scaleNum = 5;
kColumns = [10 20 35 40 50]; %linspace(10,50,scaleNum);
title_ = 'Columns';
threshold = 0.05;
dataDimPerm{1} = permute(X,[1 3 2]);

for i = 1:length(kColumns)
    [classCols{i},vectorsCols{i},affinityCols{i}] = svdClass(dataDimPerm,kColumns(i),embedded,threshold,title_,plotFlagKmeans);
end

clear dataDimPerm


% rows - apply svdClass for different scales
scaleNum = 5;
kRows = [10 20 35 40 50]; %linspace(10,50,scaleNum);
title_ = 'Rows';
threshold = 0.05;
dataDimPerm{1} = permute(X,[3 2 1]);

for i = 1:length(kRows)
    [classRows{i},vectorsRows{i},affinityRows{i}] = svdClass(dataDimPerm,kRows(i),embedded,threshold,title_,plotFlagKmeans);
end

clear dataDimPerm

if(~coupling)

    kTime = 30;
    title_ = 'Time';
    threshold = 0.05;
    dataDimPerm{1} = X;
    [classTime,vectorsTime,affinityTime] = svdClass(dataDimPerm,kTime,embedded,threshold,title_,plotFlagKmeans);
    clear dataDimPerm
    
    OrganizeData3D( dataOrig, X, affinityRows{3}, affinityCols{3}, affinityTime, rows_perm, col_perm, trial_perm, 4, 4, 4 )
    return
end

    
for j=1:length(kColumns)
        dataCoupeledTime = coupleData(permute(X,[1 3 2]),classCols{j},classRows{j});
        dataDimPerm{j} = permute(dataCoupeledTime,[1 3 2]);
end

kTime = 30;
title_ = 'Time';
threshold = 0.05;
[classTime,vectorsTime,affinityTime] = svdClass(dataDimPerm,kTime,embedded,threshold,title_,plotFlagKmeans);

clear dataDimPerm


% frames are re-organized by the embedding, columns and rows by the
% clusterring

[ err_rate,organized_data,row_order,col_order,orderTime ] = OrganizeData3D( dataOrig, X, affinityRows{3}, affinityCols{3},...
affinityTime, rows_perm, col_perm, trial_perm, 12, 12, 12 );

orderCols = getOrder(classCols{3});
orderRows = getOrder(classRows{3});
dataReordered = uint8(X(fliplr(orderRows),orderCols,orderTime));

implay(dataReordered);









if save_movies
    matrix2mov(dataOrig, 'orig.avi');
    matrix2mov(X, 'perm.avi');
    matrix2mov(dataReordered, 'PCAclust.avi');
    matrix2mov(organized_data_initmetrics, 'initMetrics.avi');
    matrix2mov(organized_data_TDtrees, 'TdPcaTrees.avi');
    matrix2mov(organized_data_BUtrees, 'BUTrees.avi');
end


