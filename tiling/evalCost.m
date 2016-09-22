function [cost, est] = evalCost(data, tiling, representation_error_fun, eps)

est = zeros(size(data));
tiles = unique(tiling(:));
for ti = 1:length(tiles)
    tileLoc = tiling == tiles(ti);
    tiled = data(tileLoc);
    est(tileLoc) = mean(tiled(:));
end
representation_error = representation_error_fun(data, est);

cost = representation_error + eps*length(unique(tiling(:)));
end
