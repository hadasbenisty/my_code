addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../SvdKmeansCluster'));

addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));


dbstop if error;
datapth = 'C:\Users\Hadas\Documents\work\datasets\circuit_nimrod\exp1\results\output';
%% init running params
ispermute = false;
dims4quest=3;
params = SetGenericDimsQuestParams(dims4quest, true);
for ind = 1:dims4quest
    params.init_aff{ind}.metric = 'euc';
    params.tree{ind}.runOnEmbdding = true;
    params.tree{ind}.treeDepth = 8;%4;
    params.tree{ind}.splitsNum = 2;
end
runningOrder = [     1 2 3 ];
n_iters=2;
params.emd{runningOrder==1}.eps=80;
params.init_aff{runningOrder==1}.knn=50;
params.tree{runningOrder == 2}.min_cluster = 6;
params.verbose=2;
params.data.over_rows = false;
params.data.to_normalize = true;%false
params.data.normalization_type = 'by_std';

%% initial analysis of the files
if ~exist('circ_params.mat','file')
    n = 1;
    pb = CmdLineProgressBar('Reading VCD files');
    
    for Y=[9 58 255]
        for M = [0 1]
            for X=0:255
                for Z=0:255
                    filename = ['X' num2str(X) '_Y' num2str(Y) '_Z' num2str(Z) '_MODE' num2str(M) '.vcd'];
                    if ~exist(fullfile(datapth, filename), 'file')
                        continue;
                    end
                    pb.print(n, 60e3);
                    circ_params.X(n) = X;
                    circ_params.Y(n) = Y;
                    circ_params.Z(n) = Z;
                    circ_params.M(n) = M;
                    n = n + 1;
                end
            end
        end
    end
    save('circ_params','circ_params');
else
    load('circ_params')
end
uniqueX = unique(circ_params.X);uniqueZ = unique(circ_params.Z);
Nparams_all = length(uniqueX);




%% first exp - X*Y+Z: mode [0 1], Y = 9, Z = oneval

NparamsX = 20;
NparamsZ= 1;
Yvals = 9;
Mvals = [0 1];
[data, wires, x, f] = loadVcdData(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, 0);
[xpermed, n1_perm, n2_perm, n3_perm] = permdata(x, ispermute);
X=data.X(f==1);
Y=data.Y(f==1);
Z=data.Z(f==1);
M=data.M(f==1);
Res = (X.*Y+Z).*(M==0)+(X.*Y+Z).*(M==1);
params.n_iters=n_iters;
[ Trees, dual_aff, initial_aff, embedding] = RunGenericDimsQuestionnaire( params, permute(xpermed,(runningOrder) ) );
for dimi = 1:3
[vecs, vals] = CalcEigs(threshold(initial_aff{dimi}, params.init_aff{dimi}.thresh)    , params.tree{dimi}.eigs_num);
initial_embedding{dimi} = vecs*vals;
end
figure;
subplot(3,2,1);plotEmbeddingWithColors(initial_embedding{runningOrder==1}(:,1:2), wires{1}.name(n1_perm), 'wires');
subplot(3,2,2);plotEmbeddingWithColors(initial_embedding{runningOrder==2}(:,1:2), n2_perm, 'Time');
subplot(3,2,3);plotEmbeddingWithColors(initial_embedding{runningOrder==3}(:,1:2), X(n3_perm), 'X');
subplot(3,2,4);plotEmbeddingWithColors(initial_embedding{runningOrder==3}(:,1:2), Z(n3_perm), 'Z');
subplot(3,2,5);plotEmbeddingWithColors(initial_embedding{runningOrder==3}(:,1:2), M(n3_perm), 'Mode');

figure;
% subplot(5,2,1);plotEmbeddingWithColors(embedding{runningOrder==1}(:,1:2), wires{1}.name(n1_perm), 'wires');
% subplot(5,2,2);plotTreeWithColors(Trees{runningOrder==1}, wires{1}.name(n1_perm));title('wires');
subplot(4,2,1);plotEmbeddingWithColors(embedding{runningOrder==2}(:,1:3), n2_perm, 'Time');
subplot(4,2,2);plotTreeWithColors(Trees{runningOrder==2}, n2_perm);title('Time');
subplot(4,2,3);plotEmbeddingWithColors(embedding{runningOrder==3}(:,1:2), X(n3_perm), 'X');
subplot(4,2,4);plotTreeWithColors(Trees{runningOrder==3},  X(n3_perm));title( 'X');
subplot(4,2,5);plotEmbeddingWithColors(embedding{runningOrder==3}(:,1:2), Z(n3_perm), 'Z');
subplot(4,2,6);plotTreeWithColors(Trees{runningOrder==3}, Z(n3_perm));title( 'Z');
subplot(4,2,7);plotEmbeddingWithColors(embedding{runningOrder==3}(:,1:2), M(n3_perm), 'Mode');
subplot(4,2,8);plotTreeWithColors(Trees{runningOrder==3}, M(n3_perm));title( 'Mode');
% suptitle(['Two Modes, Y = 9; Trees & Embedding Based on Informed Metric (' num2str(n_iters) ' Iteration)']);
figure;plotEmbeddingWithColors(embedding{runningOrder==3}(:,1:2), Res(n3_perm), 'Res');
    
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
% for n=1:3
% params.tree{n}.runOnEmbdding = false;
% end
% params.tree{runningOrder==2}.runOnEmbdding = false;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
% [~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);
% n3_orderedtree = organizeTreeForTiling(Trees{runningOrder==3}, n3_roder);


p1.eps=.2;
orderedData = xpermed(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
% orderedData = orderedData(:, :, n3_roder);

data2tile_norm =NormalizeData(orderedData, params.data);
for t2tile = [find(X(n3_roder)==16 & Z(n3_roder)==80) find(X(n3_roder)==224 & Z(n3_roder)==80)]
    
[minCurrErr, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data2tile_norm(:,:,t2tile), {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data2tile_norm(:,:,t2tile), solutionTiling);

s1 = bsxfun(@minus, data2tile_norm(:,:,t2tile).', mean(data2tile_norm(:,:,t2tile).'));
s1 = bsxfun(@rdivide, s1, std(data2tile_norm(:,:,t2tile).')+eps);

s2 = bsxfun(@minus,meanTiled.', mean(meanTiled.'));
s2 = bsxfun(@rdivide, s2, std(meanTiled.')+eps);
figure;plotTiledData(s1.', s2.',wires{1}.mdl_name, solutionTiling.isbusy,...
 ['(X,Y,Z)=(' num2str(X(n3_roder(t2tile))),',',num2str(Y(n3_roder(t2tile))),',',num2str(Z(n3_roder(t2tile))),')']);

figure;plotTiledData(data2tile_norm(:,:,t2tile), meanTiled,wires{1}.mdl_name, solutionTiling.isbusy,...
 ['(X,Y,Z)=(' num2str(X(n3_roder(t2tile))),',',num2str(Y(n3_roder(t2tile))),',',num2str(Z(n3_roder(t2tile))),')']);
% figure;
% imagesc(1:13, 1:size(orderedData, 1),orderedData(:,:,t2tile));title(['(X,Y,Z)=(' num2str(X(n1_roder(t2tile))),',',num2str(Y(n1_roder(t2tile))),',',num2str(Z(n1_roder(t2tile))),')',' Mode = ',num2str(M(n1_roder(t2tile)))]);
% set(gca, 'Ytick', 1:length(n1_roder));
% set(gca, 'YtickLabel', wiresnimrod{n1_roder});colormap gray;
% colormap gray;
end

%% Second exp - X*Y+Z: mode 0, Y = [9], X = 10 vals, Z = one val

NparamsX = 20;
NparamsZ= 1;
Yvals = [9];
Mvals = [0 1];
[data1, wires1, x1, f1] = loadVcdData(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals(1), 0);
[data2, wires2, x2, f2] = loadVcdData(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals(2), 0);

x = cat(2, x1,x2);
xn=x+0.01*randn(size(x));
X=[data1.X(f1==1) data2.X(f2==1)];
Y=[data1.Y(f1==1) data2.Y(f2==1)];
Z=[data1.Z(f1==1) data2.Z(f2==1)];
M=[data1.M(f1==1) data2.M(f2==1)];
params.n_iters=n_iters;
[ Trees, dual_aff, initial_aff, embedding] = RunGenericDimsQuestionnaire( params, permute(xn,(runningOrder) ) );
for dimi = 1:3
[vecs, vals] = CalcEigs(threshold(initial_aff{dimi}, params.init_aff{dimi}.thresh)    , params.tree{dimi}.eigs_num);
initial_embedding{dimi} = vecs*vals;
end
figure;
% subplot(3,2,1);plotEmbeddingWithColors(initial_embedding{1}(:,1:3), wires{1}.name, 'wires');
subplot(2,2,1);plotEmbeddingWithColors(initial_embedding{2}(:,1:3), 1:26, 'Time');
subplot(2,2,2);plotEmbeddingWithColors(initial_embedding{3}(:,1:2), X(1:NparamsX), 'X');
subplot(2,2,3);plotEmbeddingWithColors(initial_embedding{3}(:,1:2), Z(1:NparamsX), 'Z');
subplot(2,2,4);plotEmbeddingWithColors(initial_embedding{2}(:,1:3), [zeros(13,1); ones(13,1)], 'Mode');

figure;
subplot(4,2,1);plotEmbeddingWithColors(embedding{2}(:,1:3), 1:26, 'Time');
subplot(4,2,2);plotTreeWithColors(Trees{2}, 1:26);title('Time');
subplot(4,2,3);plotEmbeddingWithColors(embedding{3}(:,1:2), X(1:NparamsX), 'X');
subplot(4,2,4);plotTreeWithColors(Trees{3},  X(1:NparamsX));title( 'X');
subplot(4,2,5);plotEmbeddingWithColors(embedding{3}(:,1:2), Z(1:NparamsX), 'Z');
subplot(4,2,6);plotTreeWithColors(Trees{3}, Z(1:NparamsX));title( 'Z');
subplot(4,2,7);plotEmbeddingWithColors(embedding{2}(:,1:3), [zeros(13,1); ones(13,1)], 'Mode');
subplot(4,2,8);plotTreeWithColors(Trees{2}, [zeros(13,1); ones(13,1)]);title( 'Mode');
% suptitle(['Two Modes, Y = 9; Trees & Embedding Based on Informed Metric (' num2str(n_iters) ' Iteration)']);
figure;plotEmbeddingWithColors(embedding{3}(:,1:3), Res(1:NparamsX), 'Res');
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
% for n=1:3
% params.tree{n}.runOnEmbdding = false;
% end
% params.tree{runningOrder==2}.runOnEmbdding = false;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
[~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);
n3_orderedtree = organizeTreeForTiling(Trees{runningOrder==3}, n3_roder);


p1.eps=2;
orderedData = xn(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
orderedData = orderedData(:, :, n3_roder);

for t2tile = [find(X(n3_roder)==26 & Z(n3_roder)==117) find(X(n3_roder)==86 & Z(n3_roder)==117)]

    figure;
    subplot(2,2,1);
    plotTreeWithColors(n1_orderedtree, 1:length(n1_roder));
    title('wires');
    subplot(2,2,2);
    plotTreeWithColors(n2_orderedtree, 1:length(n2_roder));
    title('Time');
    subplot(2,2,3);
    imagesc(orderedData(:,:,t2tile));
    subplot(2,2,4);
[minCurrErr, solutionTiling] = reversedIterativeTilingGenericDimsleagal(orderedData(:,:,t2tile), {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(orderedData(:,:,t2tile), solutionTiling);
figure;plotTiledData(orderedData(:,:,t2tile), meanTiled,wires1{1}.mdl_name, solutionTiling.isbusy,...
 ['(X,Y,Z)=(' num2str(X(n3_roder(t2tile))),',',num2str(Y(n3_roder(t2tile))),',',num2str(Z(n3_roder(t2tile))),')']);
% figure;
% imagesc(1:13, 1:size(orderedData, 1),orderedData(:,:,t2tile));title(['(X,Y,Z)=(' num2str(X(n1_roder(t2tile))),',',num2str(Y(n1_roder(t2tile))),',',num2str(Z(n1_roder(t2tile))),')',' Mode = ',num2str(M(n1_roder(t2tile)))]);
% set(gca, 'Ytick', 1:length(n1_roder));
% set(gca, 'YtickLabel', wiresnimrod{n1_roder});colormap gray;
% colormap gray;
end

    
%%
%% Third exp - X*Y+Z: mode [0,1], Y = [9], X = 10 vals, Z = 10

NparamsX = 8;
NparamsZ= 8;
Yvals = [9];
Mvals = [0 ];
[data, wires, x, f] = loadVcdData(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, 0);

xn=x+0.01*randn(size(x));
X=data.X(f==1);
Y=data.Y(f==1);
Z=data.Z(f==1);
M=data.M(f==1);
Res = (X.*Y+Z).*(M==0)+(X.*Y+Z).*(M==1);
params.n_iters=n_iters;
[ Trees, dual_aff, initial_aff, embedding] = RunGenericDimsQuestionnaire( params, permute(xn,(runningOrder) ) );
for dimi = 1:3
[vecs, vals] = CalcEigs(threshold(initial_aff{dimi}, params.init_aff{dimi}.thresh)    , params.tree{dimi}.eigs_num);
initial_embedding{dimi} = vecs*vals;
end
figure;
% subplot(3,2,1);plotEmbeddingWithColors(initial_embedding{1}(:,1:3), wires{1}.name, 'wires');
subplot(2,2,1);plotEmbeddingWithColors(initial_embedding{2}(:,1:2), 1:13, 'Time');
subplot(2,2,2);plotEmbeddingWithColors(initial_embedding{3}(:,1:2), X, 'X');
subplot(2,2,3);plotEmbeddingWithColors(initial_embedding{3}(:,1:2), Z, 'Z');
subplot(2,2,4);plotEmbeddingWithColors(initial_embedding{3}(:,1:2), M, 'Mode');

figure;
subplot(4,2,1);plotEmbeddingWithColors(embedding{2}(:,1:2), 1:13, 'Time');
subplot(4,2,2);plotTreeWithColors(Trees{2}, 1:13);title('Time');
subplot(4,2,3);plotEmbeddingWithColors(embedding{3}(:,1:2), X, 'X');
subplot(4,2,4);plotTreeWithColors(Trees{3},  X);title( 'X');
subplot(4,2,5);plotEmbeddingWithColors(embedding{3}(:,1:2), Z, 'Z');
subplot(4,2,6);plotTreeWithColors(Trees{3}, Z);title( 'Z');
subplot(4,2,7);plotEmbeddingWithColors(embedding{3}(:,1:2), M, 'Mode');
subplot(4,2,8);plotTreeWithColors(Trees{3}, M);title( 'Mode');
% suptitle(['Two Modes, Y = 9; Trees & Embedding Based on Informed Metric (' num2str(n_iters) ' Iteration)']);
figure;plotEmbeddingWithColors(embedding{3}(:,1:3), Res, 'Res');
    

%% 
clear data, clear wires;
% load data
NparamsX = 10;
NparamsZ= 10;
Yvals = 9;
Mvals = [0,1];
[data, wires, f] = loadVcdData(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, 'mode01_y_9_xy_10vals.mat', 1);

[xpermed, n1_perm, n2_perm, n3_perm] = permdata(x);
xnpermed=xpermed+0.01*randn(size(x));




X=data.X(f==1);
Y=data.Y(f==1);
Z=data.Z(f==1);
M=data.M(f==1);
Res = (X.*Y+Z).*(M==0)+(X.*Y+Z).*(M==1);
% run quest

% initial embeddings
params.n_iters=0;
[ Trees0, ~, ~, embedding0] = RunGenericDimsQuestionnaire( params, permute(xnpermed,(runningOrder) ) );
figure;
% subplot(5,2,1);plotEmbeddingWithColors(embedding0{1}(:,1:2), wires{1}.name(n1_perm), 'wires');
% subplot(5,2,2);plotTreeWithColors(Trees0{1}, wires{1}.name(n1_perm));title('wires');
subplot(5,2,3);plotEmbeddingWithColors(embedding0{2}(:,1:2), n2_perm, 'Time');
subplot(5,2,4);plotTreeWithColors(Trees0{2}, n2_perm);title('Time');
subplot(5,2,5);plotEmbeddingWithColors(embedding0{3}(:,1:2), X(n3_perm), 'X');
subplot(5,2,6);plotTreeWithColors(Trees0{3},  X(n3_perm));title( 'X');
subplot(5,2,7);plotEmbeddingWithColors(embedding0{3}(:,1:2), Z(n3_perm), 'Z');
subplot(5,2,8);plotTreeWithColors(Trees0{3}, Z(n3_perm));title( 'Z');
subplot(5,2,9);plotEmbeddingWithColors(embedding0{3}(:,1:2), M(n3_perm), 'Mode');
subplot(5,2,10);plotTreeWithColors(Trees0{3}, M(n3_perm));title( 'Mode');
suptitle('Two Modes, Y = 9; Trees & Embedding Based on Initial Metric')

% 1 iteration
for n_iters = 1
    params.n_iters=n_iters;
    [ Trees, dual_aff, ~, embedding] = RunGenericDimsQuestionnaire( params, permute(xnpermed,(runningOrder) ) );
    figure;
%     subplot(5,2,1);plotEmbeddingWithColors(embedding{1}(:,1:2), wires{1}.name(n1_perm), 'wires');
%     subplot(5,2,2);plotTreeWithColors(Trees{1}, wires{1}.name(n1_perm));title('wires');
    subplot(5,2,3);plotEmbeddingWithColors(embedding{2}(:,1:2), n2_perm, 'Time');
    subplot(5,2,4);plotTreeWithColors(Trees{2}, n2_perm);title('Time');
    subplot(5,2,5);plotEmbeddingWithColors(embedding{3}(:,1:2), X(n3_perm), 'X');
    subplot(5,2,6);plotTreeWithColors(Trees{3},  X(n3_perm));title( 'X');
    subplot(5,2,7);plotEmbeddingWithColors(embedding{3}(:,1:2), Z(n3_perm), 'Z');
    subplot(5,2,8);plotTreeWithColors(Trees{3}, Z(n3_perm));title( 'Z');
    subplot(5,2,9);plotEmbeddingWithColors(embedding{3}(:,1:2), M(n3_perm), 'Mode');
    subplot(5,2,10);plotTreeWithColors(Trees{3}, M(n3_perm));title( 'Mode');
    suptitle(['Two Modes, Y = 9; Trees & Embedding Based on Informed Metric (' num2str(n_iters) ' Iteration)']);
    figure;plotEmbeddingWithColors(embedding{3}(:,1:2), Res(n3_perm), 'Res');
    
   
end

epsilon=[1 1];
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
% for n=1:3
% params.tree{n}.runOnEmbdding = false;
% end
% params.tree{runningOrder==2}.runOnEmbdding = false;
[~, n1_roder] = sort(Trees{1}{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{1}{runningOrder==2}{2}.clustering);
[~, n3_roder] = sort(Trees{1}{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{1}{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{1}{runningOrder==2}, n2_roder);
n3_orderedtree = organizeTreeForTiling(Trees{1}{runningOrder==3}, n3_roder);


p1.eps=.5;
orderedData = xpermed(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
orderedData = orderedData(:, :, n3_roder);

for t2tile = [find(X(n1_roder)==24 & Z(n1_roder)==81) find(X(n1_roder)==150 & Z(n1_roder)==8) find(X(n1_roder)==24 & Z(n1_roder)==190)]

%     figure;
%     subplot(2,2,1);
%     plotTreeWithColors(n1_orderedtree, 1:length(n1_roder));
%     title('wires');
%     subplot(2,2,2);
%     plotTreeWithColors(n2_orderedtree, 1:length(n2_roder));
%     title('Time');
%     subplot(2,2,3);
%     imagesc(orderedData(:,:,t2tile));
%     subplot(2,2,4);
[minCurrErr, solutionTiling] = reversedIterativeTilingGenericDimsleagal(orderedData(:,:,t2tile), {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(orderedData(:,:,t2tile), solutionTiling);
figure;plotTiledData(orderedData(:,:,t2tile), meanTiled,wiresnimrod(n1_roder), solutionTiling.isbusy,...
 ['(X,Y,Z)=(' num2str(X(n1_roder(t2tile))),',',num2str(Y(n1_roder(t2tile))),',',num2str(Z(n1_roder(t2tile))),')',' Mode = ',num2str(M(n1_roder(t2tile)))]);
% figure;
% imagesc(1:13, 1:size(orderedData, 1),orderedData(:,:,t2tile));title(['(X,Y,Z)=(' num2str(X(n1_roder(t2tile))),',',num2str(Y(n1_roder(t2tile))),',',num2str(Z(n1_roder(t2tile))),')',' Mode = ',num2str(M(n1_roder(t2tile)))]);
% set(gca, 'Ytick', 1:length(n1_roder));
% set(gca, 'YtickLabel', wiresnimrod{n1_roder});colormap gray;
% colormap gray;
end
p1.eps=.1;
p1.verbose=2;
figure;
    subplot(2,2,1);
    plotTreeWithColors(n1_orderedtree, 1:length(n1_roder));
    title('wires');
    subplot(2,2,2);
    plotTreeWithColors(n2_orderedtree, 1:length(n2_roder));
    title('Time');
    subplot(2,2,3);
    imagesc(orderedData(:,:,t2tile));
    subplot(2,2,4);
[minCurrErr, solutionTiling3d] = reversedIterativeTilingGenericDimsleagal(orderedData, {n1_orderedtree, n2_orderedtree, n3_orderedtree}, p1);

