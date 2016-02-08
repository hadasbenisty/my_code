% run the 3-D Questionnaire on flexible trees using several speakers saying 
% selected words
%% Initialization
close all;
clear all;
clc;
addpath(genpath('../Questionnaire'));
addpath(genpath('../utils'))
wavspath = '..\..\..\datasets\TIMIT\TEST';
sentenceName = 'SA1';
% set the seed for reproducibility
rng(6548964);%356454
dbstop if error;
% set parameters
run_on_mfcc = false;
isRandPerm = false;

dorandperm_col = isRandPerm;
dorandperm_row = isRandPerm;
dorandperm_height = isRandPerm;

framesize = 20e-3;
savefigs = false;
hop = 0.5;

%% Load input
selectedWords = [  1 2 ;5 6  ; 9 10     ];
dialects = dir(wavspath);
x_unpadded={};
genderNum=[];
speakerNum=[];wordsNum=[];dialectNum=[];spkrind=1;
for n = 3:length(dialects)
    speakers = dir(fullfile(wavspath, dialects(n).name));
    for m = [3 4 5 6]%:length(speakers)% [3 4 7 8] [3 8] - working nice for words separation
        [rawspeach, fs] = wavread(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WAV']));
        [wordsSamples] = readWrdFile(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WRD']), rawspeach);
        for selectedWordsind = 1:size(selectedWords, 1)
            x_unpadded{end+1}=[];
            for internalwordsinds = 1:size(selectedWords, 2)
                x_unpadded{end} = [x_unpadded{end}; wordsSamples{selectedWords(selectedWordsind, internalwordsinds)}];
            end
            if speakers(m).name(1) == 'F'
                genderNum(end+1) = 1;
            else
                genderNum(end+1) = 2;
            end
            speakerNum(end+1) = spkrind;
            
            wordsNum(end+1) = selectedWordsind;
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
    trials_perm = randperm(nT);
else
    trials_perm = 1:nT;
end
data = X(row_perm, :, :);
data = data(:, col_perm, :);
data = data(:, :, trials_perm);


%% Run Questionnaire 3D
% [row_tree, col_tree] = RunQuestionnaire(params, data(1:15, 1:20, 1));

[row_tree, col_tree, trials_tree, row_aff, col_aff, trials_aff] = RunQuestionnaire3D(params, data(:, :, :));
% close all;

%% Visualization
figspath = 'D:\workWithBoss\summaries\figs\';
coltitle =  'TimeFrames';
trialtitle = 'Trials';
rowtitle = 'Freq';
eigsnum_col = 4; eigsnum_row=4;eigsnum_trials=4;
% OrganizeData3D(X, data, row_aff, col_aff, trials_aff, row_perm, col_perm, trials_perm,  eigsnum_col, eigsnum_row, eigsnum_trials);
plotTreesAndEmbedding3D(figspath, savefigs, rowtitle, coltitle, trialtitle, ...
                                 eigsnum_row, eigsnum_col, eigsnum_trials, ...
                                 row_aff, col_aff, trials_aff, ...
                                 row_tree, col_tree, trials_tree, ...
                                 row_perm, col_perm, trials_perm, wordsNum);

% plot the trials tree
figure;
subplot(2,2,1);
plotTreeWithColors(trials_tree, genderNum(trials_perm));
title('Trials Tree By Gender');
subplot(2,2,2);
plotTreeWithColors(trials_tree,  speakerNum(trials_perm));
title('Trials Tree By Speaker');
subplot(2,2,3);
plotTreeWithColors(trials_tree,  dialectNum(trials_perm));
title('Trials Tree By Dialect');
subplot(2,2,4);
plotTreeWithColors(trials_tree,  wordsNum(trials_perm));
title('Trials Tree By Words');
if savefigs
print('D:\workWithBoss\summaries\figs\TrialsTree.pdf','-dpdf')
end
