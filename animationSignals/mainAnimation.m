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


addpath(genpath('../../PCAclustering/'));
addpath('../gen_utils/');

%% Ordering by initial metrics;
params  = SetGenericQuestParamsAnimation;

col_init_aff = feval(params.init_aff_col.initAffineFun, X, params.init_aff_col);
trial_init_aff = feval(params.init_aff_col.initAffineFun, permute(X, [2 1 3]), params.init_aff_col);
row_init_aff = feval(params.init_aff_col.initAffineFun, permute(X, [1 3 2]), params.init_aff_col);
[ err_rate2,organized_data2,row_order2,col_order2,orderTime2 ] = OrganizeData3D( dataOrig, X, trial_init_aff, col_init_aff,...
row_init_aff, rows_perm, col_perm, trial_perm, 12, 12, 12 );
organized_data_initmetrics = X(:,col_order2, :); 
organized_data_initmetrics = organized_data_initmetrics(row_order2,:, :); 
organized_data_initmetrics = organized_data_initmetrics(:,:, orderTime2);
%% run TD trees with SVD clustering
[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, X );
[ err_rate1,organized_data1,row_order1,col_order1,orderTime1 ] = OrganizeData3D( double(dataOrig), double(X), col_dual_aff, trial_dual_aff,...
row_dual_aff, rows_perm, col_perm, trial_perm, 12, 12, 12 );
% orderCols1 = getOrder(col_tree{2}.clustering);
% orderRows1 = getOrder(trial_tree{2}.clustering);
% dataReordered1 = uint8(X(fliplr(orderRows),orderCols,orderTime1));

organized_data_TDtrees = X(:,col_order1, :); 
organized_data_TDtrees = organized_data_TDtrees(row_order1,:, :); 
organized_data_TDtrees = organized_data_TDtrees(:,:, orderTime1); 


implay(uint8(dataReordered));
implay(uint8(organized_data_initmetrics));
implay(uint8(organized_data_TDtrees));
if save_movies
    matrix2mov(dataOrig, 'orig.avi');
    matrix2mov(X, 'perm.avi');
    matrix2mov(dataReordered, 'PCAclust.avi');
    matrix2mov(organized_data_initmetrics, 'initMetrics.avi');
    matrix2mov(organized_data_TDtrees, 'TdPcaTrees.avi');
end


