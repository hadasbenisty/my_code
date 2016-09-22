% run the 3-D Questionnaire on flexible trees using several speakers saying 
% selected words
%% Initialization
close all;
clear all;
clc;
addpath(genpath('../Questionnaire'));
addpath(genpath('../utils'))
wavspath = '..\..\..\datasets\RAW\TEST';
sentenceName = 'SA1';
% set the seed for reproducibility
rng(6548964);%356454
dbstop if error;
% set parameters
run_on_mfcc = false;
dorandperm_col = true;
dorandperm_row = true;
dorandperm_height = true;
framesize = 20e-3;

hop = 0.5;

%% Load input
selectedWords = [ 8 11 1 2  ];
dialects = dir(wavspath);
x_unpadded={};
genderNum=[];
speakerNum=[];wordsNum=[];dialectNum=[];spkrind=1;
for n = [3 6 8]%:length(dialects)
    speakers = dir(fullfile(wavspath, dialects(n).name));
    for m = [3 5 10]%:length(speakers)% [3 4 7 8] [3 8] - working nice for words separation
        [rawspeach, fs] = audioread(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WAV']));
        [wordsSamples] = readWrdFile(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WRD']), rawspeach);
        for selectedWordsind = selectedWords
            x_unpadded{end+1} = wordsSamples{selectedWordsind};
            if speakers(m).name(1) == 'F'
                genderNum(end+1) = 1;
            else
                genderNum(end+1) = 2;
            end
            speakerNum(end+1) = spkrind;
            
            wordsNum(end+1) = find(selectedWords == selectedWordsind);
            dialectNum(end+1) = n-2;
        end
        spkrind=spkrind+1;
    end
end
% sort according to words
[wordsNum, sortedinds] = sort(wordsNum);
speakerNum = speakerNum(sortedinds);
genderNum = genderNum(sortedinds);
dialectNum = dialectNum(sortedinds);

x_unpadded = {x_unpadded{sortedinds}};

nfft = fs * framesize;
% padd with zeros because the speach signals do share the same rate
maxsize=0;
nT = length(x_unpadded);
for T = 1:nT
    maxsize = max(maxsize, length(x_unpadded{T}));
end
xpadded = 2e-4*randn(maxsize, nT);
for T = 1:nT
    xpadded(1:length(x_unpadded{T}), T) = x_unpadded{T};
end

for T = 1:nT
    
    if run_on_mfcc
        X(:, :, T)=wav2mfcc(xpadded(:, T),fs,framesize,framesize*hop,1,1,23,22,0.1);
        eigsnum_col = 12;%20;
        eigsnum_row = 12;%75;
        
        params = getParamsForMFCC(eigsnum_col,eigsnum_row, eigsnum_trials);
    else
        % stft by Israel Cohen
        P_israel=stft(xpadded(:, T),nfft, nfft*hop, 1);
       
        X(:, :, T) = log(abs(P_israel));
        eigsnum_col = 30;%20;
        eigsnum_row = 20;%75;
        eigsnum_trials = 20;
        params = getParamsForSTFT(eigsnum_col,eigsnum_row, eigsnum_trials);
    end
    [n_rows, n_cols] = size(X(:, :, T));
end


%% perm data
if dorandperm_col
    col_perm = randperm(n_cols);
else
    col_perm = 1:n_cols;
end
if dorandperm_row
    row_perm = randperm(n_rows);
else
    row_perm = 1:n_rows;
end
if dorandperm_height
    height_perm = randperm(nT);
else
    height_perm = 1:nT;
end
data = X(row_perm, :, :);
data = data(:, col_perm, :);
data = data(:, :, height_perm);


%% Run Questionnaire 3D
% [row_tree, col_tree] = RunQuestionnaire(params, data(1:15, 1:20, 1));

[row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data(:, :, :));
% close all;

%% Visualization
% plot the trials tree
figure;
subplot(2,2,1);
plotTreeWithColors(trials_tree, genderNum(height_perm));
title('Trials Tree By Gender');
subplot(2,2,2);
plotTreeWithColors(trials_tree,  speakerNum(height_perm));
title('Trials Tree By Speaker');
subplot(2,2,3);
plotTreeWithColors(trials_tree,  dialectNum(height_perm));
title('Trials Tree By Dialect');
subplot(2,2,4);
plotTreeWithColors(trials_tree,  wordsNum(height_perm));
title('Trials Tree By Words');
print('D:\workWithBoss\summaries\figs\TrialsTree.pdf','-dpdf')

% plot the Frequency tree
figure;
freq = length(row_perm):-1:1;
plotTreeWithColors(row_tree, freq(row_perm));
title('Frequency Tree');
print('D:\workWithBoss\summaries\figs\FreqTree.pdf','-dpdf')

% plot the Time Frames tree
figure;
plotTreeWithColors(col_tree, col_perm);
title('Time Frames Tree');
print('D:\workWithBoss\summaries\figs\TimeTree.pdf','-dpdf')

 

% %% This is only for setting the thresholds       
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
% [~, i] = sort(height_perm);
% imagesc(trials_aff(i, i));colorbar;
% % keyboard;
% row_thresh = 0.9;
% col_thresh = 0.9;
% trials_thresh = 0.6;% 0.4
% row_aff = threshold(row_aff, row_thresh);
% col_aff = threshold(col_aff, col_thresh);
% trials_aff = threshold(trials_aff, trials_thresh);
%% Get Final Embedding

eigsnum_col = 4; eigsnum_row=4;eigsnum_trials=3;
% OrganizeData3D(X, data, row_aff, col_aff, trials_aff, row_perm, col_perm, height_perm,  eigsnum_col, eigsnum_row, eigsnum_trials);

[col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
embedding = col_vecs*col_vals;
PlotEmbedding( embedding, col_perm, 'Time Frames Embedding' )
print('D:\workWithBoss\summaries\figs\TimeE.pdf','-dpdf')


% draw embedding for rows
[row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
embedding = row_vecs*row_vals;
PlotEmbedding( embedding, row_perm, 'Frequencies Embedding' )
print('D:\workWithBoss\summaries\figs\FreqE.pdf','-dpdf')


[trials_vecs, trials_vals] = CalcEigs(trials_aff, eigsnum_trials);
embedding = trials_vecs*trials_vals;
color_mat = jet(2);
color_mat = color_mat(genderNum,:);
figure;
subplot(2,2,1); 
if size(embedding, 2) ==2
 scatter(embedding(:,1), embedding(:,2), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
else
    scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
end
title('Trials Embedding, Colors By Gender');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

color_mat = jet(max(speakerNum));
color_mat = color_mat(speakerNum(height_perm),:);

subplot(2,2,2); 
if size(embedding, 2) ==2
    scatter(embedding(:,1), embedding(:,2), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
else
    scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
end
title('Trials Embedding, Colors By Speaker');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
color_mat = jet(max(wordsNum));
color_mat = color_mat(wordsNum(height_perm),:);

subplot(2,2,3);
if size(embedding, 2) ==2
    scatter(embedding(:,1), embedding(:,2), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
else
    scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
end
title('Trials Embedding, Colors By Words');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

color_mat = jet(max(dialectNum));
color_mat = color_mat(dialectNum(height_perm),:);
subplot(2,2,4); 
if size(embedding, 2) ==2
    scatter(embedding(:,1), embedding(:,2), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
else
    scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
end
title('Trials Embedding, Colors By Dialect');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on
print('D:\workWithBoss\summaries\figs\EmbeddingTrials.pdf','-dpdf')
% scatter(embedding(:,1), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
% title('Trials Embedding, Colors By Words');xlabel('\psi_1'), ylabel('\psi_2'), hold on
[~, i] = sort(height_perm);
figure;imagesc(trials_aff(i,i));
colorbar;
print('D:\workWithBoss\summaries\figs\AffineTrials.pdf','-dpdf')
