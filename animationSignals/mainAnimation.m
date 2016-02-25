clc;
clear all;
close all hidden;
close all;
addpath(genpath('../../3D_Questionnaire\Questionnaire'))
addpath('../../SvdKmeansCluster');
embedded   = false;
coupling = true;
randData = false;
plotFlag = false;
plotFlagKmeans = false;
useFit = false;

load('../../SvdKmeansCluster\Animation\workspaceAnimation2');
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


