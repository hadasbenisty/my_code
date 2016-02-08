
% load wav fild
[x,fs] = wavread('arctic_a0101.wav');
% this is because fs is 16e3, for enother fs maybe other window size will
% be better
nfft=512;
% plot stft directly from matlab
spectrogram(x,nfft,'yaxis');
% get the stft from matlab
P_matlab = spectrogram(x,nfft);
% plot matlab''s stft 
figure;imagesc(log(abs(flipud(P_matlab))));
title('Matlab''s STFT');
% stft by Israel Cohen
P_israel=stft(x,nfft);
% plot Israel''s STFT
figure;imagesc(log(abs(flipud(P_israel))));
title('Israel''s STFT');
% reconstract STFT using Israel's code
x_rec=istft(Y,nfft);
% display the reconstracion error
disp(['Reconstraction Error by Israel: ' num2str(mean(abs(x(1:length(x_rec))-x_rec)))])

