addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../../SvdKmeansCluster'));

addpath(genpath('../../PCAclustering'));
addpath(genpath('../../Matlab_Utils/'));
addpath(genpath('../tiling'));
addpath(genpath('../../GraphConnectivity/'));

dbstop if error;
datapth = 'C:\Users\Hadas\Documents\work\datasets\circuit_nimrod\exp1\results_new\output';
%% init running params
ispermute = false;
dims4quest=3;
figspath='figs';
mkNewFolder(figspath);
params = SetGenericDimsQuestParams(dims4quest, true);
for ind = 1:dims4quest
    params.init_aff{ind}.metric = 'euc';
    params.tree{ind}.runOnEmbdding = true;
    params.tree{ind}.treeDepth = 8;%4;
    params.tree{ind}.splitsNum = 2;
end
overwrite = false;
runningOrder = [    2  1  3 ];
n_iters=3;
params.emd{runningOrder==1}.eps=80;
params.init_aff{runningOrder==1}.eps = 1e3;
params.emd{runningOrder==1}.eps = 1e3;
params.init_aff{runningOrder==1}.knn=50;
params.tree{runningOrder == 2}.min_cluster = 6;
params.tree{runningOrder == 1}.eigs_num=30;
params.verbose=1;
params.data.over_rows = true;
params.data.to_normalize = true;%false
params.data.normalization_type = 'by_std';

%% initial analysis of the files
circ_params = initial_circ_files_analysis('circ_params.mat', datapth);
uniqueX = unique(circ_params.X);uniqueZ = unique(circ_params.Z);
Nparams_all = length(uniqueX);

%% load data
%% A
NparamsX = 40;
NparamsZ= 1;
Yvals = 58;
Mvals = [0 1];
isbin2dec = false;
mainLabels = {'u_adder'    'u_adderOut'    'u_cnt'    'u_fsm'    'u_muxA'    'u_muxB'    'u_muxOut' 'u_muxXorZ'};

[data, wires, samples, f] = loadVcdData_moduled(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, overwrite, isbin2dec);
[Labels, numLabels] = getModuleLabels(wires{1}.module, mainLabels);
X=data.X(f==1);
Y=data.Y(f==1);
Z=data.Z(f==1);
M=data.M(f==1);
Res = (X.*Y+Z).*(M==0)+(X.*Y+Z).*(M==1);
mode0indic = M==0;
mode1indic = M==1;
%% first exp - X = 40 values, mode = 0, Y = 9 or 58, Z = oneval
x = samples(:,:,mode0indic);
x_cumsum = cumsum(x,2);
params.n_iters=0;
figure;imagesc(x(:,:,1));colormap gray; xlabel('Time[samples]');
ylabel('Wires'); mysave(gcf, fullfile(figspath, 'exp1_samples'));

[ Trees0, ~, initial_aff0] = RunGenericDimsQuestionnaire( params, permute(x_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp1_initial', {'Module', 'Time', 'X'}, [3 3  3], ...
    initial_aff0(runningOrder), -inf*[1  1 1], ...
    Trees0(runningOrder), {[] []  []}, {numLabels, 1:62, X(mode0indic)});
params.n_iters=n_iters;
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(x_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp1_final', {'Module', 'Time', 'X'}, [3 3 3 ], ...
    dual_aff(runningOrder), -inf*[1 1 1 ], ...
    Trees(runningOrder), {[] [] [] }, {numLabels, 1:62, X(mode0indic)});




% tiling
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
% [~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);


p1.eps=12000;
orderedData = x(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
data4tiling = mean(orderedData,3);
data4tilingLog = log(threshold(data4tiling+2*abs(min(data4tiling(:))),0)+eps);

if p1.verbose>1
    figure;
    subplot(2,2,3);plotTreeWithColors(Trees{runningOrder==1}, numLabels(n1_roder));title('wires');view(-90,90);
    subplot(2,2,2);plotTreeWithColors(Trees{runningOrder==2}, 1:62);title('Time');
    subplot(2,2,4);
end
[~, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data4tilingLog, {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data4tiling, solutionTiling);
f=figure;
plotTiledData(data4tilingLog, meanTiled,Labels(n1_roder), solutionTiling.isbusy,'Exp. 1' );
mysave(gcf, fullfile(figspath, 'exp1_tiled'));
mysave(f, fullfile(figspath, 'exp1_blueTiles'));

figure;
imagesc(1:13, 1:size(orderedData, 1),data4tilingLog);title('Exp. 1');
set(gca, 'Ytick', 1:length(n1_roder));
set(gca, 'YtickLabel', Labels);colormap gray;
colormap gray;
mysave(gcf, fullfile(figspath, 'exp1_data'));



%% second exp - X*Y+Z: mode [0 1], Y = 9 or 58, Z = oneval concatenating the mode
ca;
x = samples;
x_cumsum = cumsum(x,2);
figure;imagesc(x(:,:,1));colormap gray; xlabel('Time[samples]');
ylabel('Wires'); mysave(gcf, fullfile(figspath, 'exp2_samples'));

params.n_iters=0;
[ Trees0, ~, initial_aff0] = RunGenericDimsQuestionnaire( params, permute(x_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp2_initial', {'Module', 'Time', 'X', 'Mode'}, [3 3 3 3], ...
    initial_aff0([runningOrder find(runningOrder==3)]), -inf*[1 1 1 1], ...
    Trees0([runningOrder find(runningOrder==3)]), {[] [] [] []}, {numLabels, 1:62, X, M});
params.n_iters=2;
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(x_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp2_final', {'Module', 'Time', 'X', 'Mode'}, [3 3 3 3], ...
    dual_aff([runningOrder find(runningOrder==3)]), -inf*[1 1 1 1], ...
    Trees([runningOrder find(runningOrder==3)]), {[] [] [] []}, {numLabels, 1:62, X, M});
% tiling
p1.representation_error_fun = @squaredErr;
p1.verbose = 2;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
% [~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);


p1.eps=12000;
orderedData = x(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
data4tiling = mean(orderedData,3);
data4tilingLog = log(threshold(data4tiling+2*abs(min(data4tiling(:))),0)+eps);

if p1.verbose>1
    figure;
    subplot(2,2,3);plotTreeWithColors(Trees{runningOrder==1}, numLabels(n1_roder));title('wires');view(-90,90);
    subplot(2,2,2);plotTreeWithColors(Trees{runningOrder==2}, 1:62);title('Time');
    subplot(2,2,4);
end
[~, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data4tilingLog, {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data4tiling, solutionTiling);
f=figure;
plotTiledData(data4tilingLog, meanTiled,Labels(n1_roder), solutionTiling.isbusy,'Exp. 2' );
mysave(gcf, fullfile(figspath, 'exp2_tiled'));
mysave(f, fullfile(figspath, 'exp2_blueTiles'));

figure;
imagesc(1:13, 1:size(orderedData, 1),data4tilingLog);title('Exp. 2');
set(gca, 'Ytick', 1:length(n1_roder));
set(gca, 'YtickLabel', Labels);colormap gray;
colormap gray;
mysave(gcf, fullfile(figspath, 'exp2_data'));

%% Exp 3 

x_twomodes_concat = [x(:,:,mode0indic) x(:,:,mode1indic)];
x_twomodes_concat_cumsum = cumsum(x_twomodes_concat,2);
figure;imagesc(x_twomodes_concat(:,:,1));colormap gray; xlabel('Time[samples]');
ylabel('Wires'); mysave(gcf, fullfile(figspath, 'exp3_samples'));

X_twomodes_concat = X(mode0indic);
Y_twomodes_concat = Y(mode0indic);
Z_twomodes_concat = Z(mode0indic);
M_twomodes_concat = [zeros(62, 1); ones(62, 1)];
[Labels, numLabels] = getModuleLabels(wires{1}.module, mainLabels);

params.n_iters=0;
[ Trees0, ~, initial_aff0] = RunGenericDimsQuestionnaire( params, permute(x_twomodes_concat_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp3_initial', {'Module', 'Time', 'X', 'Mode'}, [3 3 3 3], ...
    initial_aff0([runningOrder find(runningOrder==2)]), -inf*[1 1 1 1], ...
    Trees0([runningOrder find(runningOrder==2)]), {[] [] [] []}, {numLabels, 1:124, X_twomodes_concat, M_twomodes_concat});


params.n_iters=2;
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(x_twomodes_concat_cumsum,(runningOrder) ) );
PlotTreesAndEmbedding(true, figspath, 'exp3_final', {'Module', 'Time', 'X', 'Mode'}, [3 3 3 3], ...
    dual_aff([runningOrder find(runningOrder==2)]), -inf*[1 1 1 1], ...
    Trees([runningOrder find(runningOrder==2)]), {[] [] [] []}, {numLabels, 1:124, X_twomodes_concat, M_twomodes_concat});

% tiling
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
[~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);


p1.eps=12e3;
orderedData = x_twomodes_concat(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
data4tiling = mean(orderedData,3);
data4tilingLog = log(threshold(data4tiling+2*abs(min(data4tiling(:))),0)+eps);

if p1.verbose>1
    figure;
    subplot(2,2,3);plotTreeWithColors(Trees{runningOrder==1}, numLabels(n1_roder));title('wires');view(-90,90);
    subplot(2,2,2);plotTreeWithColors(Trees{runningOrder==2}, 1:124);title('Time');
    subplot(2,2,4);
end
[minCurrErr, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data4tilingLog, {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data4tiling, solutionTiling);
f=figure;
plotTiledData(data4tilingLog, meanTiled,Labels(n1_roder), solutionTiling.isbusy,'Exp. 3' );
mysave(gcf, fullfile(figspath, 'exp3_tiled'));
mysave(f, fullfile(figspath, 'exp3_blueTiles'));

figure;
imagesc(1:13, 1:size(orderedData, 1),data4tilingLog);title('Exp. 3');
set(gca, 'Ytick', 1:length(n1_roder));
set(gca, 'YtickLabel', Labels);colormap gray;
colormap gray;
mysave(gcf, fullfile(figspath, 'exp3_data'));

%% B
NparamsX = 4;
NparamsZ= 40;
Yvals = 58;
Mvals = [0 1];
isbin2dec = false;

[data, wires, samples, f] = loadVcdData_moduled(datapth, uniqueX, uniqueZ, Nparams_all, NparamsX, Yvals, NparamsZ, Mvals, overwrite, isbin2dec);
[Labels, numLabels] = getModuleLabels(wires{1}.module, mainLabels);
X=data.X(f==1);
Y=data.Y(f==1);
Z=data.Z(f==1);
M=data.M(f==1);
Res = (X.*Y+Z).*(M==0)+(X.*Y+Z).*(M==1);
uniqex = unique(X);
for u = 1:length(uniqex)
    xindic(:,u) = (X == uniqex(u) & M == 1);
end
%% Exp. 4 - X,Z = 40 values, mode = 0, Y = 58
x_4modes_concat = [samples(:,:,xindic(:,1)) samples(:,:,xindic(:,2)) samples(:,:,xindic(:,3)) samples(:,:,xindic(:,4))];
x_4modes_concat_cumsum = cumsum(x_twomodes_concat,2);
figure;imagesc(x_4modes_concat(:,:,1));colormap gray; xlabel('Time[samples]');
ylabel('Wires'); mysave(gcf, fullfile(figspath, 'exp4_samples'));

X_4modes_concat = kron(uniqex,ones(1,62));
Y_4modes_concat = Y(xindic(:,1));
Z_4modes_concat = Z(xindic(:,1));
M_4modes_concat = M(xindic(:,1));

params.n_iters=0;
[ Trees0, ~, initial_aff0] = RunGenericDimsQuestionnaire( params, permute(x_4modes_concat,(runningOrder) ) );
PlotTreesAndEmbedding(false, figspath, 'exp4_initial', {'Module', 'Time',  'Z', 'M','X'}, [3 3  3 3 3], ...
    initial_aff0([runningOrder find(runningOrder==3) find(runningOrder==2)]), -inf*[1  1 1 1 1], ...
    Trees0([runningOrder find(runningOrder==3) find(runningOrder==2)]), {[] []  [] [] []}, {numLabels, 1:248, ...
    Z_4modes_concat, M_4modes_concat, X_4modes_concat});
params.n_iters=2;
[ Trees, dual_aff] = RunGenericDimsQuestionnaire( params, permute(x_4modes_concat,(runningOrder) ) );
PlotTreesAndEmbedding(false, figspath, 'exp4_final', {'Module', 'Time',  'Z', 'M','X'}, [3 3  3 3 3], ...
    dual_aff([runningOrder find(runningOrder==3) find(runningOrder==2)]), -inf*[1  1 1 1 1], ...
    Trees([runningOrder find(runningOrder==3) find(runningOrder==2)]), {[] []  [] [] []}, {numLabels, 1:248, ...
    Z_4modes_concat, M_4modes_concat, X_4modes_concat});
[~, neurons_order_all] = sort(Trees{runningOrder==1}{2}.clustering);



% tiling
p1.representation_error_fun = @squaredErr;
p1.verbose = 0;
[~, n1_roder] = sort(Trees{runningOrder==1}{2}.clustering);
n2_roder = getTimeOrderFromTree(Trees{runningOrder==2}{2}.clustering);
% [~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Trees{runningOrder==1}, n1_roder);
n2_orderedtree = organizeTreeForTiling(Trees{runningOrder==2}, n2_roder);


p1.eps=5000;
orderedData = x_twomodes_concat(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
data4tiling = orderedData(:,:,1);
data4tilingLog = log(threshold(data4tiling+2*abs(min(data4tiling(:))),0)+eps);

if p1.verbose>1
    figure;
    subplot(2,2,3);plotTreeWithColors(Trees{runningOrder==1}, numLabels(n1_roder));title('wires');view(-90,90);
    subplot(2,2,2);plotTreeWithColors(Trees{runningOrder==2}, 1:62);title('Time');
    subplot(2,2,4);
end
[~, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data4tilingLog, {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data4tiling, solutionTiling);
f=figure;
plotTiledData(data4tilingLog, meanTiled,Labels(n1_roder), solutionTiling.isbusy,'Exp. 1' );
mysave(gcf, fullfile(figspath, 'exp1_tiled'));
mysave(f, fullfile(figspath, 'exp1_blueTiles'));

figure;
imagesc(1:13, 1:size(orderedData, 1),data4tilingLog);title('Exp. 1');
set(gca, 'Ytick', 1:length(n1_roder));
set(gca, 'YtickLabel', Labels);colormap gray;
colormap gray;
mysave(gcf, fullfile(figspath, 'exp1_data'));

% 
% 
% params.connectivity.winLen = 1;
% params.connectivity.overlap = .5;
% params.connectivity.epsParamD = 1;
% params.connectivity.epsParamVis = 1;
% params.connectivity.distTypeAlltimes = 3;
% params.connectivity.distTypeTimeSample = 4;
% params.connectivity.meanTrials = false;
% 
% [ DallTimes_time, D_time, Vis_time, VisNorm_time, Vis_Dalltimes_time ] = ...
%     GraphConnectivity( permute(x_twomodes_concat_cumsum, [1 2 3]),  params.connectivity);
% 
% [ DallTimes_permuted, D_neurons_permuted, Vis_neurons_permuted, VisNorm_neurons_permuted, Vis_Dalltimes_neurons_permuted ] = ...
%     GraphConnectivity( permute(x_twomodes_concat_cumsum, [2 1 3]),  params.connectivity);
% 
% [ DallTimes_trials, D_trials, Vis_trials, VisNorm_trials, Vis_Dalltimes_trials ] = ...
%     GraphConnectivity( permute(x_twomodes_concat_cumsum, [1 3 2]),  params.connectivity);
% 
% 
% 
% DallTimes_neurons = DallTimes_permuted(neurons_order_all, neurons_order_all);
% D_neurons = D_neurons_permuted(:,:,neurons_order_all );
% D_neurons = D_neurons(:,neurons_order_all,: );
% Vis_neurons = Vis_neurons_permuted(:,:, neurons_order_all);
% Vis_neurons = Vis_neurons(:,neurons_order_all, :);
% Vis_Dalltimes_neurons = Vis_Dalltimes_neurons_permuted(neurons_order_all, neurons_order_all);
% 

