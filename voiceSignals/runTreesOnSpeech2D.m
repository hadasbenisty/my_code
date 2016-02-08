% run the 3-D Questionnaire on flexible trees using several utterances of the same sentence
%% Initialization
close all;
clear all;
clc;
addpath(genpath('../Questionnaire'));
addpath(genpath('../STFT'))
addpath(genpath('voicebox'));
wavspath = 'D:\workWithBoss\datasets\TIMIT\TEST';
sentenceName = 'SA1';
% set the seed for reproducibility
rng(73631);
dbstop if error;
% set parameters
run_on_mfcc = false;
dorandperm_col = true;
dorandperm_row = true;
framesize = 20e-3;

hop = 0.3;

% Load input
selectedWords = [1 2 ];
sentenceind=1;
dialects = dir(wavspath);
x_unpadded={};
genderNum=[];
speakerNum=[];wordsNum=[];voicingNum=[];vadNum=[];
for n = 3:3%length(dialects)
    speakers = dir(fullfile(wavspath, dialects(n).name));
    for m = 3:length(speakers)
        [rawspeach, fs] = wavread(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WAV']));
        [wordsSamples] = readWrdFile(fullfile(wavspath, dialects(n).name, speakers(m).name, [sentenceName '.WRD']), rawspeach);
        
        x=[];wordsNumcurr=[];voicingcurr=[];vadcurr=[];
        for selectedWordsind = selectedWords
            voicingcurr=[voicingcurr v_ppmvu(wordsSamples{selectedWordsind},fs,'au').'];
            vadcurr = [vadcurr vadsohn(wordsSamples{selectedWordsind},fs,'b').'];
            x = [x; wordsSamples{selectedWordsind}];
            wordsNumcurr = [wordsNumcurr ones(1, length(wordsSamples{selectedWordsind}))*selectedWordsind];
        end
        wordsNum{end+1} = wordsNumcurr;
        voicingNum{end+1} = voicingcurr;
        vadNum{end+1} = vadcurr;
        x_unpadded{end+1} = x;
        if speakers(m).name(1) == 'F'
            genderNum(end+1) = 1;
        else
            genderNum(end+1) = 2;
        end
        speakerNum(end+1) = m-2;
        
    end
end

nfft = fs * framesize;
X = [];wordsNumsFinal=[];voicingNumFinal=[];vadNumFinal=[];
for T = 1:length(x_unpadded)
    if run_on_mfcc
        error('not supported yet');
        X(:, :, T)=wav2mfcc(x_unpadded{T},fs,framesize,framesize*hop,1,1,23,22,0.1);
        eigsnum_col = 12;%20;
        eigsnum_row = 12;%75;
        
        params = getParamsForMFCC(eigsnum_col,eigsnum_row);
    else
        % stft by Israel Cohen
        P_israel=stft(x_unpadded{T},nfft, nfft*hop, 1);
        wordsNumsampled = wordsNum{T}(1:nfft*hop:end); 
        wordsNumsampled=wordsNumsampled(1:size(P_israel, 2));
        voicingNumSampled= voicingNum{T}(1:nfft*hop:end); 
        voicingNumSampled=voicingNumSampled(1:size(P_israel, 2));
        vadNumSampled= vadNum{T}(1:nfft*hop:end); 
        vadNumSampled=vadNumSampled(1:size(P_israel, 2));
        
        
        X = [X  log(abs(P_israel))];
        wordsNumsFinal = [wordsNumsFinal wordsNumsampled];
        voicingNumFinal = [voicingNumFinal voicingNumSampled];
        vadNumFinal = [vadNumFinal vadNumSampled];
        eigsnum_col = 30;%20;
        eigsnum_row = 20;%75;
        
        params = getParamsForSTFT(eigsnum_col,eigsnum_row);
    end
    [n_rows, n_cols] = size(X);
end


% perm data
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

data = X(row_perm, :, :);
data = data(:, col_perm, :);


%% Run Questionnaire 3D
% [row_tree, col_tree] = RunQuestionnaire(params, data(1:15, 1:20, 1));

[row_tree, col_tree] = RunQuestionnaire(params, data);
close all;

%% Visualization
% plot the col. tree
% figure, treeplot(nodes(col_tree),'.')
% 
% mannercolors = zeros(size(nodes(col_tree)));
% mannercolors(length(mannercolors)-n_cols+1:length(mannercolors)) = manner_num;
% [x_layout,y_layout] = treelayout(nodes(col_tree));
% k=1;
% color_mat = jet(length(unique(manner_num)));
% color_mat = color_mat(manner_num+1,:);
% hold on;
% for i=length(x_layout)-n_cols+1:length(x_layout)
%     text(x_layout(i),y_layout(i)-0.02,manner_str{k});
%     d=plot(x_layout(i),y_layout(i),'.');
%     set(d,'MarkerEdgeColor',color_mat(k,:),'MarkerFaceColor',color_mat(k,:),'MarkerSize',10);
%     hold on;
%     k=k+1;
% end
% Get Final Embedding
 
row_aff = CalcEmdAff(data.', col_tree, params.row_emd, params.init_aff_row.on_rows);

col_aff = CalcEmdAff(data, row_tree, params.col_emd, ~params.init_aff_col.on_rows);

        
% figure;subplot(2,2,1);
% hist(row_aff(:),100);
% subplot(2,2,3);
% hist(col_aff(:),100);
% subplot(2,2,2);
% imagesc(row_aff);colorbar;
% subplot(2,2,4);
% imagesc(col_aff);colorbar;

% keyboard;
row_thresh = 0.5;
col_thresh = 0.5;
% fairly high thresholds on the affinities (data-dependent), which are aimed
% to collapse the diffusion map into a curve.
row_aff = threshold(row_aff, row_thresh);
col_aff = threshold(col_aff, col_thresh);

eigsnum_col = 30; eigsnum_row=20;
% OrganizeData((X), data, row_aff, col_aff, row_perm, col_perm, eigsnum_col, eigsnum_row);
% draw embedding for cols
[col_vecs, col_vals] = CalcEigs(col_aff, eigsnum_col);
embedding = col_vecs*col_vals;

color_mat = jet(size(embedding, 1));
color_mat = color_mat(col_perm,:);

figure, scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Time Frames Embedding');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

color_mat = jet(max(wordsNumsFinal));
color_mat = color_mat(wordsNumsFinal,:);
figure;scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Time Frames Embedding, Color By Words');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

voicingNumFinal=voicingNumFinal/max(voicingNumFinal);
color_mat = jet(256);
color_mat = color_mat(1+round(voicingNumFinal*255),:);
figure;scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Time Frames Embedding, Color By Voicing');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

vadNumFinal=vadNumFinal/max(vadNumFinal);
color_mat = jet(256);
color_mat = color_mat(1+round(vadNumFinal*255),:);

figure;scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Time Frames Embedding, Color By VAD');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on


[~,i]=sort(col_perm);
figure;subplot(6,1,1);
imagesc(X);
subplot(6,1,2);
plot(embedding(i,1));
subplot(6,1,3);
plot(embedding(i,2));
subplot(6,1,4);
plot(embedding(i,3));
subplot(6,1,5);
plot(voicingNumFinal);
subplot(6,1,6);
plot(vadNumFinal);
% draw embedding for rows
[row_vecs, row_vals] = CalcEigs(row_aff, eigsnum_row);
embedding = row_vecs*row_vals;

color_mat = jet(length(embedding(:,1)));
color_mat = color_mat(row_perm,:);

figure; scatter3(embedding(:,1), embedding(:,2), embedding(:,3), 30, color_mat, 'filled'), %'MarkerEdgeColor', 'k'
title('Frequencies Embedding');xlabel('\psi_1'), ylabel('\psi_2'), zlabel('\psi_3'), hold on

