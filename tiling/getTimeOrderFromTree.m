function time_order = getTimeOrderFromTree(clustering)

clusters = unique(clustering);
for ci = 1:length(clusters)
   inds(ci) = mean(find(clustering == clusters(ci))); 
end
[~, iinds] = sort(inds, 'ascend');
newclustering = zeros(size(clustering));
for ci = 1:length(clusters)
   newclustering(clustering == clusters(iinds(ci))) = ci; 
end

[~, time_order] = sort(newclustering);
