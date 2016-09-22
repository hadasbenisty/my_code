function plotTiledData(orderedData, meanTiled, NeuronsLabels, tiling, ttl, t)

if nargin == 5
    t = 1:size(orderedData, 2);
end
delt = 0.5;

% imagesc(t, 1:size(orderedData, 1), meanTiled);colormap gray;
% hold all;
% placeTiles(tiling, delt, t)
% if ~isempty(NeuronsLabels)
% set(gca, 'Ytick', 1:length(NeuronsLabels));
% set(gca, 'YtickLabel', NeuronsLabels);
% end
% figure;
% subplot(1,2,1);
imagesc(t, 1:size(orderedData, 1), orderedData);title(ttl);
hold all;
if ~isempty(NeuronsLabels)
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels);
end
% delt = (t(2)-t(1));
hold all;
placeTiles(tiling, delt, t)
colormap gray;
figure;
% % subplot(1,2,2);
imagesc(t, 1:size(orderedData, 1), meanTiled);colormap gray;
hold all;
placeTiles(tiling, delt, t)

title(['Tiled Data of ' ttl]);
if ~isempty(NeuronsLabels)
set(gca, 'Ytick', 1:length(NeuronsLabels));
set(gca, 'YtickLabel', NeuronsLabels);
end


end

function placeTiles(tiling, delt, t)
delt_t = (t(2)-t(1))*delt;
tiles = unique(tiling(:));
for ci = 1:length(tiles)
   [ind_i, ind_j] = find(tiling == tiles(ci)) ;
   if length(ind_i) == 1 && length(ind_j) == 1
       plot(ind_j, ind_i, 'bx');
   end
   line((t(min(ind_j))-delt_t)*[1 1],[min(ind_i)-delt max(ind_i)+delt ],'Color','b', 'LineWidth',2);
   line((t(max(ind_j))+delt_t)*[1 1],[min(ind_i)-delt max(ind_i)+delt  ],'Color','b', 'LineWidth',2);  
   line([t(min(ind_j))-delt_t t(max(ind_j))+delt_t ], (min(ind_i)-delt)*[1 1],'Color','b', 'LineWidth',2); 
   line([t(min(ind_j))-delt_t t(max(ind_j))+delt_t ],(max(ind_i)+delt)*[1 1],'Color','b', 'LineWidth',2); 
end
end
