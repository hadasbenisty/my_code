function plotByClustering(tree, data, labels, newLabels, ttl)

meanMat=[];allMat=[];
folders = unique(tree.clustering);
for ci = 1:length(folders)
    inds2folder = find(tree.clustering == folders(ci));
    currLabels = labels(inds2folder);
    loc = findStrLocInCellArray(newLabels, currLabels);
    meanMat(ci, :) = mean(mean(permute(data(loc(loc~=0), :, :), [2 3 1]),2), 3);
    allMat{ci} = permute(mean(permute(data(loc(loc~=0), :, :), [2 3 1]),2), [3 1 2]);
end
figure;
imagesc(meanMat);
set(gca, 'Ytick', 1:length(folders));title(ttl);
figure;
for k=1:length(allMat)
   subplot(length(allMat), 1, k);
   imagesc(allMat{k});
   set(gca, 'XtickLabel',{});
   set(gca, 'YtickLabel',{});
   ylabel(num2str(k));
end
