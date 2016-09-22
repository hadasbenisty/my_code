close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../../3D_Questionnaire/Questionnaire'));
addpath(genpath('../../3D_Questionnaire/utils'));
addpath(genpath('../tiling'));


figspath1 = fullfile('tilingDemoOutput');


rng(654164);

%% Init params

eigsnum_col = 2;
eigsnum_row = 3;
eigsnum_trial= 2;
row_alpha = .2;
row_beta = 0;
col_alpha = .2;
col_beta = 0;
trial_alpha = .2;
trial_beta = 0;
params  = SetQuest3DParams(eigsnum_col, eigsnum_row, eigsnum_trial, row_alpha, row_beta, col_alpha, col_beta, trial_alpha, trial_beta );
params.init_aff_row_metric = 'euc';
params.init_aff_col_metric = 'euc';
%% Data

data = [1 1  1   1  2  2  3  3  3  3;...
    1 1  1   1  2  2  3  3  3  3;...
    4 4  8   8  11 11 11 11 80 80;...
    4 4  8   8  12 12 12  12 80 80;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 18 18 40 40 40 40;...
    6 14 16  17 70 70 70 70 90 90;...
    6 14 16  17 95 95 95 95 90 90];
data = data + randn(size(data))*0.3;
figure;
imagesc(data);colorbar;
title('Data');

%% Run Qu. 2D
[row_tree, col_tree, row_dual_aff, col_dual_aff] = RunQuestionnaire(params, data);

%% Visualization
[vecs_col, vals_col] = CalcEigs(col_dual_aff, 4);
[vecs_row, vals_row] = CalcEigs(row_dual_aff, 4);
subplot(2,1,1);plotEmbeddingWithColors(vecs_col*vals_col, 1:size(vecs_col,1), 'Col')
subplot(2,1,2);plotEmbeddingWithColors(vecs_row*vals_row, 1:size(vecs_row,1), 'Row')
[ row_order ] = OrganizeDiffusion3DbyOneDim( data, vecs_row*vals_row );
[ col_order ] = OrganizeDiffusion3DbyOneDim( data.', vecs_col*vals_col );

% [~, row_order] = sort(row_tree{2}.clustering);
% [~, col_order] = sort(col_tree{2}.clustering);
orderedData = data(row_order, :);
orderedData = orderedData(:, col_order);


% prepare trees for recursion
for treeLevel = 1:length(row_tree)
    row_orderedtree{treeLevel} = row_tree{treeLevel};
    row_orderedtree{treeLevel}.clustering = row_tree{treeLevel}.clustering( row_order);
end
for treeLevel = 1:length(col_tree)
    col_orderedtree{treeLevel} = col_tree{treeLevel};
    col_orderedtree{treeLevel}.clustering = col_tree{treeLevel}.clustering( col_order);
end
figure;
subplot(2,2,1);plotTreeWithColors(row_tree, 1:length(row_dual_aff));    title('Time Tree');
subplot(2,2,2);plotTreeWithColors(col_tree, 1:length(col_dual_aff));    title('Freq Tree');
subplot(2,2,3);plotTreeWithColors(row_orderedtree, 1:length(row_dual_aff));    title('Time Tree');
subplot(2,2,4);plotTreeWithColors(col_orderedtree, 1:length(col_dual_aff));    title('Freq Tree');

figure;
subplot(2,2,1);
imagesc(data);title('Data');
subplot(2,2,2);
imagesc(orderedData);title('Ordered Data');
subplot(2,2,3);
plotTreeWithColors(row_orderedtree, 1:size(data,1));
title('Row Tree');

subplot(2,2,4);
plotTreeWithColors(col_orderedtree, 1:size(data,2))
title('Col Tree');


ind2data = 1;
minErr = Inf;
tiling.isbusy = zeros(size(orderedData, 1), size(orderedData, 2), 1);
tiling.isLeader = zeros(size(orderedData, 1), size(orderedData, 2), 1);
solutionTiling = [];
figure;
ind2data = 1;minVal=4;
figure;
maxVal = 1e3;
p1.vol = getVolsFromTrees(col_orderedtree, row_orderedtree, minVal, maxVal);
p1.vol = setdiff(p1.vol, size(data).'*[1:100]);
p1.vol = sort(p1.vol, 'ascend');
p1.vol = sort(p1.vol, 'descend');

p1.beta_col = 0;
p1.beta_row = 0;
p1.verbose = 1;
%             p1.err_fun = @evalTilingErr;
p1.err_fun = @evalTilingErrTauMeas;
tau = 0.1:0.1:1.9;
for ti = 1:length(tau)
p1.tau = tau(ti);
normData = orderedData/sqrt(sum(sum(orderedData.^2)));

[minCurrErr(ti), currSolutionTiling{ti}] = loopTiling2DEfficient(normData, row_orderedtree, col_orderedtree, p1);
[err(ti), meanTiled{ti}, taumeas(ti), b(ti)] = evalTilingErrTauMeas(normData, currSolutionTiling{ti}, p1);
figure;plotTiledData(normData, meanTiled{ti}, [], currSolutionTiling{ti}.isbusy, num2str(p1.tau))

end
%% Evaluate the coefficients and their barriers

d = zeros(numel(normData));
for n1 = 1:numel(normData)
    disp(n1);
    for n2 = 1:numel(normData)
        
        [I1,J1] = ind2sub(size(normData),n1);
        [I2,J2] = ind2sub(size(normData),n2);
        if I1 == J1 || I2 == J2
            continue;
        end
        for ki=1:length(col_orderedtree)
            if row_orderedtree{ki}.clustering(I1)==row_orderedtree{ki}.clustering(I2)
                Ni = row_orderedtree{ki}.folder_sizes( row_orderedtree{ki}.clustering(I2));
                break;
            end
        end
        for kj=1:length(col_orderedtree)
            if col_orderedtree{kj}.clustering(J1)==col_orderedtree{kj}.clustering(J2)
                Nj = col_orderedtree{kj}.folder_sizes( col_orderedtree{kj}.clustering(J2));
                break;
            end
        end
        d(n1, n2) = Ni * Nj;
    end
end

distmat = squareform(pdist(normData(:), 'cityblock'));
d(d==inf) = 0;

for ti = 1:length(tau)
tiles = unique(currSolutionTiling{ti}.isbusy(:));
N = length(tiles);
phi = zeros([size(data) N]);
for n = 1:N
    [ind1, inds2] = find(currSolutionTiling{ti}.isbusy == tiles(n));
   phi(ind1, inds2,n) =  1;
   A{ti}(n) = (sum(sum(phi(:,:,n))));  
   phi(:,:,n) = phi(:,:,n)/sqrt(sum(sum(phi(:,:,n).^2))); 
   c{ti}(n) = sum(sum(phi(:,:,n).*normData));
end

al = 0.1:0.1:1;
for alind = 1:length(al)
    
    C = (distmat./(d.^al(alind)));
C(isinf(C)) = 0;
CH(ti, alind) = max(C(:));
B(ti, alind) = CH(ti, alind)/N^(1/tau(ti)-0.5)*(sum(A{ti}.^tau(ti)))^(1/tau(ti));
end
% figure;plot(B)
% hold all;
% plot(c, 'k');title('Coefficients (black) and C_H*A_k^{\alpha+0.5}');

end