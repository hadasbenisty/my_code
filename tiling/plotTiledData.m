function plotTiledData(orderedData, meanTiled, NeuronsLabels, tiling, ttl)



subplot(1,3,1);imagesc(orderedData);title(ttl);
if ~isempty(NeuronsLabels)
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels);
end
for ci = 1:length(unique(tiling(:)))
   [ind_i, ind_j] = find(tiling == ci) ;
   line(min(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k') 
   line(max(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k')  
   line([min(ind_j) max(ind_j) ], min(ind_i)*[1 1],'Color','k') 
   line([min(ind_j) max(ind_j) ],max(ind_i)*[1 1],'Color','k') 
end
subplot(1,3,2);imagesc(tiling);
title(['Tiling of ' ttl]);
if ~isempty(NeuronsLabels)
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels);
end
% for ci = 1:length(unique(tiling(:)))
%    [ind_i, ind_j] = find(tiling == ci) ;
%    line(min(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k') 
%    line(max(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k')  
%    line([min(ind_j) max(ind_j) ], min(ind_i)*[1 1],'Color','k') 
%    line([min(ind_j) max(ind_j) ],max(ind_i)*[1 1],'Color','k') 
% end
subplot(1,3,3);imagesc(meanTiled);
title(['Tiled Data of ' ttl]);
if ~isempty(NeuronsLabels)
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels);
end

% for ci = 1:length(unique(tiling(:)))
%    [ind_i, ind_j] = find(tiling == ci) ;
%    line(min(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k') 
%    line(max(ind_j)*[1 1],[min(ind_i) max(ind_i) ],'Color','k')  
%    line([min(ind_j) max(ind_j) ], min(ind_i)*[1 1],'Color','k') 
%    line([min(ind_j) max(ind_j) ],max(ind_i)*[1 1],'Color','k') 
% end

