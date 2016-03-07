clc;
clear all;
close all hidden;
close all;

% addpath(genpath('..\3DQuest\3DQuest\Questionnaire'));
embedded   = false;
dataType = 'Animation'; % 'probField' 'Animation' 'chirp'
coupling = true;
randData = false;
plotFlag = false;
plotFlagKmeans = false;
useFit = false;
sizeRows = 300;
sizeCols = 400;
maxTime  = 200;

i = linspace(0,2*pi,sizeRows);
j = linspace(0,2*pi,sizeCols);
t = linspace(0.2,1.2,maxTime);
[J,I,T] = meshgrid(j,i,t);

if (randData)
    switch dataType
        case 'Animation'
            load('Data3D/DataAnimation')
            % implay(dataOrig)
            [sizeRows,sizeCols,maxTime] = size(dataOrig);
        case 'probField'
            probField = 0.5*(1+sin((2.*I+2.*J+5*T+3*I.*J.*T)/2));
            % probField = getWave(sizeRows,sizeCols,maxTime);
            clear I J T
            dataOrig    = bsxfun(@gt,probField,rand(sizeRows,sizeCols,maxTime));
            dataOrig    = double(dataOrig)*2-1;
        case 'chirp'
            dataOrig = getChirp;
            [sizeRows,sizeCols,maxTime] = size(dataOrig);

    end
   
    shuffleRows = randperm(sizeRows);
    shuffleCols = randperm(sizeCols);
    shuffleTime = randperm(maxTime);
    data = double(dataOrig(shuffleRows,shuffleCols,shuffleTime));

else
     %load('Data3D\workspaceAnimation3Compact0');
     %[sizeRows,sizeCols,maxTime] = size(data);
     
%     filePrefix = 'D30/8_17_14_46-80'; %0.55
%     filePrefix = 'D30/8_17_14_1-45'; %0.55
    filePrefix = '../../datasets\biomed\D30/8_12_14_1-40'; % 0.5
%     filePrefix = 'D30/8_15_13_1-35'; %0.5

    
        
%     load(strcat(filePrefix,'matrix.mat'));
%     load(strcat(filePrefix,'sensors_dat_titles.mat'));
%     load('D30/commonIndex')
%     matrix_8_15_13 = matrix(index8_15_13,:,:);
%     filePrefix = 'D30/8_12_14_1-40';
%     load(strcat(filePrefix,'matrix.mat'));
%     load(strcat(filePrefix,'sensors_dat_titles.mat'));
%     matrix_8_12_14 = matrix(index8_12_14,:,:);
%     
%     matrix = zeros([size(matrix_8_15_13,1) size(matrix_8_15_13,2) size(matrix_8_15_13,3)+size(matrix_8_12_14,3)]);
%     matrix(:,:,1:size(matrix_8_15_13,3)) = matrix_8_15_13;
%     matrix(:,:,size(matrix_8_15_13,3)+1:end) = matrix_8_12_14;
    [matrix, expLabel, NeuronsLabels] = loadNeuronsData('../../../datasets\biomed\D30', {'8_12_14_1-40'}, 360);
    [sizeRows,sizeCols,maxTime] = size(matrix);
    shuffleRows = 1:sizeRows;
    shuffleCols = 1:sizeCols;
    shuffleTime = 1:maxTime;
    data = matrix;
    dataOrig = data;
end     


% compute principal components and classes of Rows
% title_ = 'Rows';
% threshold = 0.05;
% hirarchNum = 5;
% % k = linspace(10,50,hirarchNum);
% k=10;
% dataDimPerm{1} =permute(data,[3 2 1]);
% for i = 1:length(k);
%     [classRows{i},vectorsRows{i},affinityRows{i}] = svdClass(dataDimPerm,k(i),embedded,threshold,title_,plotFlagKmeans);
% end


iter = 1;
figure;
for k_t = 5%2:7
for k_r = 8%5:10
    
% compute principal components and classes of Time
%k_t = 5;
title_ = 'Time';
threshold = 0.05;
dataDimPerm{1} = data;
[classTime,vectorsTime,affinityTime] = svdClass(dataDimPerm,k_t,embedded,threshold,title_,plotFlagKmeans);
clear dataDimPerm

% compute principal components and classes of cloumns
hirarchNum = 5;
k_c = linspace(10,50,hirarchNum);
%k_c=35;
title_ = 'Columns';
threshold = 0.7;
dataDimPerm{1} = permute(data,[1 3 2]);
for i = 1:length(k_c)
    [classCols{i},vectorsCols{i},affinityCols{i}] = svdClass(dataDimPerm,k_c(i),embedded,threshold,title_,plotFlagKmeans);
end
clear dataDimPerm

if(~coupling)
    % compute principal components and classes of cloumns
%     k_t = 5;
%     title_ = 'Time';
%     threshold = 0.05;
%     dataDimPerm{1} = data;
%     [classTime,vectorsTime,affinityTime,eigValsTime] = svdClass(dataDimPerm,k_t,embedded,threshold,title_,plotFlagKmeans);
%     showEmbedding(vectorsTime(:,2:4),eigValsTime,shuffleTime,title_);

    % compute principal components and classes of Rows
    k_r = 10;
    title_ = 'Rows';
    threshold = 0.05;
    dataDimPerm{1} = permute(data,[3 2 1]);
    [classRows,vectorsRows,affinityRows] = svdClass(dataDimPerm,k_r,embedded,threshold,title_,plotFlagKmeans);
    clear dataDimPerm
    
    OrganizeData3D( dataOrig, data, affinityRows, affinityCols{3}, affinityTime, shuffleRows, shuffleCols, shuffleTime, 4, 4, 4 )
    return
end

for i=1:1
    
%     dataCoupeledCols = coupleData(data,classTime,classRows);
% 
%     % compute principal components and classes of cloumns
%     k = 35;
%     title_ = 'Columns';
%     threshold = 0.05;
%     dataDimPerm{1} = permute(dataCoupeledCols,[1 3 2]);
%     [classCols,vectorsCols,affinityCols,eigValsCols] = svdClass(dataDimPerm,k,embedded,threshold,title_,plotFlagKmeans);
    % showEmbedding(vectorsCols(:,2:4),eigValsCols,shuffleCols,title_);
 
    
    
%     for j=1:length(k)
%             dataCoupeledTime = coupleData(permute(data,[1 3 2]),classCols{j},classRows{j});
%             dataDimPerm{j} = permute(dataCoupeledTime,[1 3 2]);
%     end
%    
%     % compute principal components and classes of Time
%     k = [10 20 30 40];
%     title_ = 'Time';
%     for j = 1:length(k)
%     [classTime{j},vectorsTime{j},affinityTime{j}] = svdClass(dataDimPerm,k(j),embedded,threshold,title_,plotFlagKmeans);
%     end
    
    

%     clear dataDimPerm
%     k = [10 20 30 40];
%     title_ = 'Time';
%     for j=1:length(k)
%         dataCoupeledCols = coupleData(data,classTime{j},classRows{j});
%         dataDimPerm{j} = permute(dataCoupeledCols,[1 3 2]);
%     end
%     k=35;
%     [classCols,vectorsCols,affinityCols] = svdClass(dataDimPerm,k,embedded,threshold,title_,plotFlagKmeans);

    for j=1:length(k_c)
            dataCoupeledTime = coupleData(permute(data,[2 1 3]),classTime,classCols{j});
            dataDimPerm{j} = permute(dataCoupeledTime,[1 3 2]);
    end
    % compute principal components and classes of Rows
    %k_r = 10;
    title_ = 'Rows';
    threshold = 0.55;
    [classRows,vectorsRows,affinityRows] = svdClass(dataDimPerm,k_r,embedded,threshold,title_,plotFlagKmeans);
    
    if(i==1)
        plotFlagKmeans = false;
    end
end


threshold = 0.4;%[0.4 0.45 0.5 0.55 0.6 0.65];%0.05:0.05:0.95;

for i = 1:length(threshold)
    
    embedded = true;
    [classRows,vectorsRows,affinityRows] = svdClass(dataDimPerm,k_r,embedded,threshold(i),title_,plotFlagKmeans);
    
    classNum = length(unique(classRows));
    colors = jet(classNum);
    color_mat = zeros(length(classRows),3);
    
    for j =1:classNum
        classIndex = find(classRows == j);
        classSize = length(classIndex);
        color_mat(classRows == j,:) = repmat(colors(j,:),[classSize 1]);

    end    
    
    [row_vecs, row_vals] = CalcEigs(affinityRows, 4);
    embedding = row_vecs*row_vals;
    subplot(6, 6, iter)
    scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled');
    title(sprintf('Neurons Embedding thresh=%.1f k_r=%.1f k_t=%.1f',threshold(i),k_r,k_t));
    xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
    xlim([-0.1 0.3])
    xlim([-0.05 0.1])
    view([-115 10])
    iter = iter + 1;
    
end
clear dataDimPerm

end
end

% filePrefix = strcat('D30/ClassResults2',filePrefix(4:end)); %0.5

% save(strcat(filePrefix,'affinity'),'affinityRows','classRows');
% save(strcat(filePrefix,'affinityCols'),'affinityCols','classCols');

% [ err_rate,row_order,col_order,trial_order ] = OrganizeData3D( dataOrig, data, affinityRows, affinityCols{1},...
%  affinityTime, shuffleRows, shuffleCols, shuffleTime, 4, 4, 4 );


