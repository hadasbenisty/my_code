function [currCost, tilingSol] = reversedIterativeTiling2Dleagal(data, row_tree, col_tree, params)

tilingSol.isbusy = ones(size(data));
currCost = evalCost(data, tilingSol.isbusy, params.eps);
MAX_ITERS = 100;

for iters = 1:MAX_ITERS
    disp(iters);
    tiles = unique(tilingSol.isbusy(:));
    splitting = zeros(length(tiles), 1);
    for ti = 1:length(tiles)
%         [inds_r, inds_c] = find(tilingSol.isbusy == tiles(ti));
        [row_split, col_split] = getNextSplitFromTree(row_tree, col_tree, tilingSol.isbusy, tiles(ti));
        if length(unique(row_split)) == 1
            cost_r = inf;
        else
            tilingSol_r.isbusy = splitTile(tilingSol.isbusy, tiles(ti), row_split, ones(size(col_split)));
            cost_r = evalCost(data, tilingSol_r.isbusy, params.eps);
        end
        if length(unique(col_split)) == 1
            cost_c = inf;
        else
            tilingSol_c.isbusy = splitTile(tilingSol.isbusy, tiles(ti), ones(size(row_split)), col_split);
            cost_c = evalCost(data, tilingSol_c.isbusy, params.eps);
        end
        [~, mini] = min([currCost, cost_r, cost_c]);
        switch mini
            case 3
                currCost = cost_c;
                tilingSol = tilingSol_c;
                meanTiled = getTiledData(data, tilingSol);
                imagesc(meanTiled);
                splitting(ti) = 1;
                drawnow;
            case 2
                currCost = cost_r;
                tilingSol = tilingSol_r;
                meanTiled = getTiledData(data, tilingSol);
                imagesc(meanTiled);
                drawnow;
                splitting(ti) = 1;
            case 1
                % do nothing
        end
disp(currCost);
    end
    if all(splitting == 0)
        break;
    end
end
disp(iters);
end
function [row_split, col_split] = getNextSplitFromTree(row_tree, col_tree, isbusy, tile)

[inds_r, inds_c] = find(isbusy == tile);
len_r = length(unique(inds_r));
len_c = length(unique(inds_c));
for tree_level_r = 1:length(row_tree)
    if length(unique(row_tree{tree_level_r}.clustering(inds_r(1):inds_r(1) + len_r -1))) == 1
        break;
    end
end
for tree_level_c = 1:length(col_tree)
    if length(unique(col_tree{tree_level_c}.clustering(inds_c(1):inds_c(1) + len_c -1))) == 1
        break;
    end
end
if tree_level_r == 1
    row_split = 1;
else
    row_split = row_tree{tree_level_r-1}.clustering(inds_r(1):inds_r(1) + len_r -1);
end
if tree_level_c == 1
    col_split = 1;
else
    col_split = col_tree{tree_level_c-1}.clustering(inds_c(1):inds_c(1) + len_c -1);
end

end

function newIsbusy = splitTile(isbusy, tile, row_split, col_split)
newIsbusy = isbusy;
[inds_r, inds_c] = find(isbusy == tile);
tile_inds_r = min(inds_r):max(inds_r);
tile_inds_c = min(inds_c):max(inds_c);
clusters_r = unique(row_split);
clusters_c = unique(col_split);
for ci_r = 1:length(clusters_r)
    for ci_c = 1:length(clusters_c)
        
        newIsbusy(tile_inds_r(find(row_split == clusters_r(ci_r))), tile_inds_c(find(col_split == clusters_c(ci_c)))) = max(newIsbusy(:)) + 1;
    end
end

end
