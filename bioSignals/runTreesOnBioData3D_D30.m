close all;
clear all;
clc;
dbstop if error;
addpath(genpath('../Questionnaire'));
rng(73631);
%% Init params
dorandperm_col = false;
dorandperm_row = false;
dorandperm_trials = false;
eigsnum_col = 100;
eigsnum_row = 100;
eigsnum_trials = 20;
params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trials);
%% Load Data

filefirst = '..\..\..\datasets\biomed\D30\8_17_14_1-45matrix.mat';
filesecond = '..\..\..\datasets\biomed\D30\8_17_14_46-80matrix.mat';

load(filefirst);
nt = 360;
nT1 = size(matrix, 2) / nt;
for T = 1:nT1
   X(:, :,  T) = matrix(:, 1 + nt*(T-1):nt * T);   
end

load(filesecond);
nT2 = size(matrix, 2) / nt;
for T = 1:nT2
   X(:, :,  T+nT1) = matrix(:, 1 + nt*(T-1):nt * T);   
end
beforeafterLabels = [ones(nT1, 1); ones(nT2, 1)*2];
[nr, nt, nT] = size(X);
% perm data
if dorandperm_col
    col_perm = randperm(nt);
else
    col_perm = 1:nt;
end
if dorandperm_row
    row_perm = randperm(nr);
else
    row_perm = 1:nr;
end
if dorandperm_trials
    trial_perm = randperm(nT);
else
    trial_perm = 1:nT;
end
data = X(row_perm, :, :);
data = data(:, col_perm, :);
data = data(:, :, trial_perm);
%% Run Qu.

[row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data);


%% Visualization

% plot the trials tree
figure;
subplot(1,2,1);
plotTreeWithColors(trials_tree, (trial_perm));
title('Trials Tree By T');

subplot(1,2,2);
plotTreeWithColors(trials_tree, beforeafterLabels(trial_perm));
title('Trials Tree By Virus');



init_trial_affin = CalcInitAff3D(permute(data, [1 3 2]), params.init_aff_row);
eigsnum_trials=12;
% OrganizeData((X), data, row_aff, col_aff, row_perm, col_perm, eigsnum_col, eigsnum_row);
% draw embedding for cols
[col_vecs, col_vals] = CalcEigs(init_trial_affin, eigsnum_trials);
embedding = col_vecs*col_vals;
color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(trial_perm,:);

figure; scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Initial Trials Embedding, Colors By Trial Index');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

trials_thresh = 0.2;% 0.4
% fairly high thresholds on the affinities (data-dependent), which are aimed
% to collapse the diffusion map into a curve.
trials_aff1 = threshold(trials_aff, trials_thresh);

[col_vecs, col_vals] = CalcEigs(trials_aff1, eigsnum_trials);
embedding = col_vecs*col_vals;
color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(trial_perm,:);

figure; scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Final Trials Embedding, Colors By Trial Index');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
xlim([-.4 0.3])
ylim([-.3 0])
zlim([-0.3 1]);

color_mat = jet(2);
color_mat = color_mat(beforeafterLabels,:);
figure; scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Final Trials Embedding, Colors By Trial Index');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
xlim([-.4 0.3])
ylim([-.3 0])
zlim([-0.3 1]);
% subplot(2,2,2);
% plotTreeWithColors(trials_tree,  speakerNum(trial_perm));
% title('Trials Tree By Speaker');
% subplot(2,2,3);
% plotTreeWithColors(trials_tree,  wordsNum(trial_perm));
% title('Trials Tree By Dialect');
% subplot(2,2,4);
% plotTreeWithColors(trials_tree,  wordsNum(trial_perm));
% title('Trials Tree By Words');
% print('TrialsTree.pdf','-dpdf')
% 
% % plot the row tree
% figure;
% freq = length(row_perm):-1:1;
% plotTreeWithColors(row_tree, freq(row_perm));
% title('Frequency Tree');
% print('FreqTree.pdf','-dpdf')
% 
% % plot the col tree
% figure;
% plotTreeWithColors(col_tree, col_perm);
% title('Time Frames Tree');
% print('TimeTree.pdf','-dpdf')
% 
%  

       
% figure;subplot(3,2,1);
% hist(row_aff(:),100);
% subplot(3,2,3);
% hist(col_aff(:),100);
% subplot(3,2,5);
% hist(trials_aff(:),100);
% subplot(3,2,2);
% [~, i] = sort(row_perm);
% imagesc(row_aff(i, i));colorbar;
% subplot(3,2,4);
% [~, i] = sort(col_perm);
% imagesc(col_aff(i, i));colorbar;
% subplot(3,2,6);
% [~, i] = sort(trial_perm);
% imagesc(trials_aff(i, i));colorbar;
% keyboard;
row_thresh = 0.0;
col_thresh = 0.0;
trials_thresh = 0.2;% 0.4
% fairly high thresholds on the affinities (data-dependent), which are aimed
% to collapse the diffusion map into a curve.
row_aff1 = threshold(row_aff, row_thresh);
col_aff1 = threshold(col_aff, col_thresh);
trials_aff1 = threshold(trials_aff, trials_thresh);
% Get Final Embedding

eigsnum_col = 12; eigsnum_row=12;eigsnum_trials=12;
% OrganizeData((X), data, row_aff, col_aff, row_perm, col_perm, eigsnum_col, eigsnum_row);
% draw embedding for cols
[col_vecs, col_vals] = CalcEigs(col_aff1, eigsnum_col);
embedding = col_vecs*col_vals;

color_mat = jet(size(embedding, 1));
color_mat = color_mat(col_perm,:);

figure,subplot(1,2,1); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Embedding by t');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
tonetime = 42;
toneLabel = [ones(tonetime, 1); 2*ones(size(col_vecs, 1)-tonetime, 1)];
color_mat = jet(2);
color_mat = color_mat(toneLabel,:);
subplot(1,2,2); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Embedding by t, Colors by tone');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on


% draw embedding for rows
[row_vecs, row_vals] = CalcEigs(row_aff1, eigsnum_row);
embedding = row_vecs*row_vals;

color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(row_perm,:);

figure; scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Embedding by r');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on


[trials_vecs, trials_vals] = CalcEigs(trials_aff1, eigsnum_trials);
embedding = trials_vecs*trials_vals;

color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(trial_perm,:);

% subplot(3,2,4); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
% title('Trials Embedding, Colors By Trial Index');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
% xlim([-1 1])
% ylim([-1 1])
% 
% color_mat = jet(2);
% color_mat = color_mat(beforeafterLabels,:);
% 
% subplot(3,2,5); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
% title('Trials Embedding, Colors By Virus Injection');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
% % xlim([-.3 0])
% % ylim([-.3 0])