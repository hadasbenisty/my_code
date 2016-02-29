function plotByClustering(meanMat, allMat, ttl, timeFrames)

if nargin == 3
    timeFrames = 1:size(meanMat, 2);
end
figure;
subplot(1,2,1);
imagesc(meanMat(:, timeFrames));
set(gca, 'Ytick', 1:size(meanMat,2));title(sprintf('%s\nMean Image', ttl));
for k=1:length(allMat)
   subplot(length(allMat), 2, 2*k);
   imagesc(allMat{k}(:,timeFrames));
   set(gca, 'XtickLabel',{});
   set(gca, 'YtickLabel',{});
   colorbar;
   ylabel(num2str(k));
   if k == 1
       title('Ordered Data');
   end
end
% figure;
% for k=1:size(meanMat,1)
%    subplot(length(allMat), 1, k);
%    imagesc(meanMat(k, :));
%    set(gca, 'XtickLabel',{});
%    set(gca, 'YtickLabel',{});
%    ylabel(num2str(k));
% end
