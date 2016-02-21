function plotByClustering(meanMat, allMat, ttl)


figure;
subplot(1,2,1);
imagesc(meanMat);
set(gca, 'Ytick', 1:size(meanMat,2));title(sprintf('%s\nMean Image', ttl));
for k=1:length(allMat)
   subplot(length(allMat), 2, 2*k);
   imagesc(allMat{k});
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
