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

load('../../../datasets/Animation\workspaceAnimation2');
[sizeRows,sizeCols,maxTime] = size(data);

% columns - apply svdClass for different scales
scaleNum = 5;
kColumns = [10 20 35 40 50]; %linspace(10,50,scaleNum);
title_ = 'Columns';
threshold = 0.05;
dataDimPerm{1} = permute(data,[1 3 2]);

for i = 1:length(kColumns)
    [classCols{i},vectorsCols{i},affinityCols{i}] = svdClass(dataDimPerm,kColumns(i),embedded,threshold,title_,plotFlagKmeans);
end

clear dataDimPerm


% rows - apply svdClass for different scales
scaleNum = 5;
kRows = [10 20 35 40 50]; %linspace(10,50,scaleNum);
title_ = 'Rows';
threshold = 0.05;
dataDimPerm{1} = permute(data,[3 2 1]);

for i = 1:length(kRows)
    [classRows{i},vectorsRows{i},affinityRows{i}] = svdClass(dataDimPerm,kRows(i),embedded,threshold,title_,plotFlagKmeans);
end

clear dataDimPerm

if(~coupling)

    kTime = 30;
    title_ = 'Time';
    threshold = 0.05;
    dataDimPerm{1} = data;
    [classTime,vectorsTime,affinityTime] = svdClass(dataDimPerm,kTime,embedded,threshold,title_,plotFlagKmeans);
    clear dataDimPerm
    
    OrganizeData3D( dataOrig, data, affinityRows{3}, affinityCols{3}, affinityTime, shuffleRows, shuffleCols, shuffleTime, 4, 4, 4 )
    return
end

    
for j=1:length(kColumns)
        dataCoupeledTime = coupleData(permute(data,[1 3 2]),classCols{j},classRows{j});
        dataDimPerm{j} = permute(dataCoupeledTime,[1 3 2]);
end

kTime = 30;
title_ = 'Time';
threshold = 0.05;
[classTime,vectorsTime,affinityTime] = svdClass(dataDimPerm,kTime,embedded,threshold,title_,plotFlagKmeans);

clear dataDimPerm


% frames are re-organized by the embedding, columns and rows by the
% clusterring

[ err_rate,organized_data,row_order,col_order,orderTime ] = OrganizeData3D( dataOrig, data, affinityRows{3}, affinityCols{3},...
affinityTime, shuffleRows, shuffleCols, shuffleTime, 4, 4, 4 );

orderCols = getOrder(classCols{3});
orderRows = getOrder(classRows{3});
dataReordered = uint8(data(fliplr(orderRows),orderCols,orderTime));

implay(dataReordered);

%% run TD trees with SVD clustering
addpath(genpath('../../PCAclustering/'));
addpath('../gen_utils/');
params  = SetGenericQuestParamsAnimation;
% para
% ms.row_emd.beta=.5;
% params.col_emd.beta=.5;
% params.trial_emd.beta=.5;
col_init_aff = feval(params.init_aff_col.initAffineFun, data, params.init_aff_col);
trial_init_aff = feval(params.init_aff_col.initAffineFun, permute(data, [2 1 3]), params.init_aff_col);
row_init_aff = feval(params.init_aff_col.initAffineFun, permute(data, [1 3 2]), params.init_aff_col);
[ err_rate2,organized_data2,row_order2,col_order2,orderTime2 ] = OrganizeData3D( dataOrig, data, trial_init_aff, col_init_aff,...
row_init_aff, shuffleRows, shuffleCols, shuffleTime, 4, 3, 4 );



[ col_tree, trial_tree, row_tree,  col_dual_aff, trial_dual_aff, row_dual_aff] = RunGenericQuestionnaire3D( params, data );
[ err_rate1,organized_data1,row_order1,col_order1,orderTime1 ] = OrganizeData3D( dataOrig, data, col_dual_aff, trial_dual_aff,...
row_dual_aff, shuffleRows, shuffleCols, shuffleTime, 4, 3, 4 );
orderCols1 = getOrder(col_tree{2}.clustering);
orderRows1 = getOrder(trial_tree{2}.clustering);
dataReordered1 = uint8(data(fliplr(orderRows),orderCols,orderTime1));

implay(uint8(organized_data1));
implay(uint8(organized_data));

organized_data111 = data(:,col_order1, :); 

organized_data111 = organized_data111(row_order1,:, :); 
organized_data111 = organized_data111(:,:, orderTime1); 

[row_vecs, row_vals] = CalcEigs(col_dual_aff, 3);
v = [ones(1,104) 100*ones(1,232)];
figure;PlotEmbedding(row_vecs*row_vals, v(shuffleRows), 'Row Embedding');



movDataOrig = matrix2mov(dataOrig);movie2avi(movDataOrig, 'orig.avi');
movPerm = matrix2mov(data);movie2avi(movPerm, 'perm.avi');
movSVMclust = matrix2mov(dataReordered);movie2avi(movSVMclust, 'PCAclust.avi');
movTDtrees = matrix2mov(organized_data1);movie2avi(movTDtrees, 'TdPcaTrees.avi');


MOV1 = aviread('TdPcaTrees.avi');


vidObj = VideoWriter('TdPcaTrees.avi');
    open(vidObj);
 
    for k = 1:length(MOV1)
       
       writeVideo(vidObj,MOV1(k));
    end
  
    % Close the file.
    close(vidObj);