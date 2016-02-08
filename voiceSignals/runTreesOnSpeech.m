% run the 2-D Questionnaire on flexible trees using one sentence
%% Initialization
close all;
clear all;
clc;
addpath(genpath('../Questionnaire'));
addpath(genpath('../STFT'))
wavspath = 'D:\workWithBoss\datasets\TIMIT\TEST';
% set the seed for reproducibility
rng(73631);
dbstop if error;
% set parameters
run_on_mfcc = false;
dorandperm_col = true;
dorandperm_row = true;
framesize = 20e-3;
nfft=round(framesize*fs);
hop = 0.3;

% Load input
[x,fs] = wavread(fullfile(wavspath, '\DR1\FAKS0\SA2.WAV'));
phones = readPhnFile(fullfile(wavspath, '\DR1\FAKS0\SA2.PHN'), nfft, hop);


if run_on_mfcc
    orig_data=wav2mfcc(x,fs,framesize,framesize*hop,1,1,23,22,0.1);
    eigsnum_col = 12;%20;
    eigsnum_row = 12;%75;
    
    params = getParamsForMFCC(eigsnum_col,eigsnum_row);
else
    % stft by Israel Cohen
    P_israel=stft(x,nfft, nfft*hop, 1);
    % reconstract STFT using Israel's code
    x_rec=istft(P_israel,nfft);
    orig_data = log(abs(P_israel));
    eigsnum_col = 12;%20;
    eigsnum_row = 12;%75;
    
    params = getParamsForSTFT(eigsnum_col,eigsnum_row);
end
% perm data
[n_rows, n_cols] = size(orig_data);
if dorandperm_col
    col_perm = randperm(n_cols);
else
    col_perm=1:n_cols;
end
if dorandperm_row
    row_perm = randperm(n_rows);
else
    row_perm=1:n_rows;
end
data = orig_data(row_perm,:);
data = data(:,col_perm);
% get metadata for later visualization
phones = {phones{1:n_cols}};
[phones_str, phones_num]  = standardTIMITphones61to39(phones);
[manner_str, manner_num] = mapPhones2Manner(phones);
[vad_str, vad_num] = mapPhones2Vad(phones);
[voiced_str, voiced_num] = mapPhones2voicedUnvoices(phones);

%% Run Questionnaire
[row_tree, col_tree] = RunQuestionnaire(params, data);
close all;

%% Visualization
% plot the col. tree
figure, treeplot(nodes(col_tree),'.')

mannercolors = zeros(size(nodes(col_tree)));
mannercolors(length(mannercolors)-n_cols+1:length(mannercolors)) = manner_num;
[x_layout,y_layout] = treelayout(nodes(col_tree));
k=1;
color_mat = jet(length(unique(manner_num)));
color_mat = color_mat(manner_num+1,:);
hold on;
for i=length(x_layout)-n_cols+1:length(x_layout)
    text(x_layout(i),y_layout(i)-0.02,manner_str{k});
    d=plot(x_layout(i),y_layout(i),'.');
    set(d,'MarkerEdgeColor',color_mat(k,:),'MarkerFaceColor',color_mat(k,:),'MarkerSize',10);
    hold on;
    k=k+1;
end
% Get Final Embedding
row_aff = CalcEmdAff(data.', col_tree, params.row_emd);
col_aff = CalcEmdAff(data, row_tree, params.col_emd);
figure;subplot(2,2,1);
hist(row_aff(:),100);
subplot(2,2,2);
hist(col_aff(:),100);
subplot(2,2,3);
imagesc(row_aff);colorbar;
subplot(2,2,4);
imagesc(col_aff);colorbar;
% keyboard;
row_thresh = 0.9;
col_thresh = 0.9;
% fairly high thresholds on the affinities (data-dependent), which are aimed
% to collapse the diffusion map into a curve.
row_aff = threshold(row_aff, row_thresh);
col_aff = threshold(col_aff, col_thresh);
eigsnum_col = 12; eigsnum_row=12;
OrganizeData((orig_data), data, row_aff, col_aff, row_perm, col_perm, eigsnum_col, eigsnum_row);
% draw embedding for cols
[col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
embedding = col_vecs*col_vals;

color_mat = jet(length(unique(manner_num)));
color_mat = color_mat(manner_num+1,:);

figure,subplot(3,2,1); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Col. Embedding, Colors By Manner');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

color_mat = jet(41);
color_mat = color_mat(phones_num,:);

subplot(3,2,2); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Col. Embedding, Colors By Phones');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

color_mat = jet(2);
color_mat = color_mat(vad_num+1,:);

subplot(3,2,3); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Col. Embedding, Colors By Vad');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
color_mat = jet(2);
color_mat = color_mat(voiced_num+1,:);

subplot(3,2,4); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Col. Embedding, Colors By Voiced');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

% draw embedding for rows
[row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
embedding = row_vecs*row_vals;

color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(row_perm,:);

subplot(3,2,5); scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Row Embedding');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

