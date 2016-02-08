function [minErr, tiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data, tiling, minErr)

MAX_S = numel(data);
if ind2data > MAX_S ||  sum(tiling.isbusy(:) == 0) == 0
    currErr = evalTilingErr(data, tiling);
    if currErr < minErr
        minErr =  currErr;
        disp(['Volume = ' num2str(vol) ' err = ' num2str(minErr)]);
        disp(tiling.isbusy);    
        disp('---------------------');
        return;
    end
end
[row_i,col_i] = ind2sub(size(data),ind2data);
if tiling.isbusy(row_i, col_i) % check if it's busy and shift to the next one
    [minErr, tiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data+1, tiling, minErr);
    return
end

tilesList = getTiles(vol);
for ti = 1:length(tilesList)
    % cheking if this tile fits
    isTileFits = checkTile(data, row_tree, col_tree, ind2data, tilesList(ti));
    if isTileFits
        % marking this tile 
        tiling = markTile(tiling, tilesList(ti), ind2data);
        % see if we don't already get an error higher than the prev. so
        % there is a point to continue
        currErr = evalTilingErr(data, tiling);
        if currErr < minErr
            [minErr, tiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data+1, tiling, minErr);
        end
        % remove this tile and continue to examine the other tiles in the
        % list
        tiling = unmarkTile(tiling, tilesList(ti), ind2data);
    end
end
end

function tiling = markTile(tiling, tile, ind2data)
       
%          /* put tile on ind2data */
%          mark ind2data as leader of tile;
[row_i,col_i] = ind2sub(size(tiling.isbusy),ind2data);
tiling.isLeader(row_i,col_i) = true;
tiling.isbusy(row_i:row_i+tile.row-1, col_i:col_i+tile.col-1) = max(tiling.isbusy(:))+1;
imagesc(tiling.isbusy);
title(['Tiling Vol = ' num2str(tile.row*tile.col)]);
drawnow;

%          mark all neighbors M[sx .. sx+tx-1][sy .. sy+ty-1] = busy;
end

function tiling = unmarkTile(tiling, tile, ind2data)
[row_i,col_i] = ind2sub(size(tiling.isbusy),ind2data);
tiling.isLeader(row_i,col_i) = false;
tiling.isbusy(row_i:row_i+tile.row-1, col_i:col_i+tile.col-1) = 0;
imagesc(tiling.isbusy);title(['Tiling Vol = ' num2str(tile.row*tile.col)]);
drawnow;

end

function isTileFits = checkTile(data, row_tree, col_tree, ind2data, tile)

% first assuming false
isTileFits = false;

[size_row, size_col] = size(data);
[row_i,col_i] = ind2sub(size(data),ind2data);

% checking if the tile is too long in each dim
if row_i + tile.row - 1 > size_row || col_i + tile.col - 1 > size_col
    return;
end

isTileFits = checkTileByTree(row_tree, tile.row, row_i);
if ~isTileFits
    return;
end
isTileFits = checkTileByTree(col_tree, tile.col, col_i);
end

function isTileFitsByTree = checkTileByTree(tree, tile_row_size, loc)

isTileFitsByTree = false;

if tile_row_size == 1
    isTileFitsByTree = true;
    return;
end

indsOfTile = loc:loc+tile_row_size-1;
for treeLevel = 1:length(tree)
    folders = tree{treeLevel}.clustering(indsOfTile);
    %check if we reached a single folder
    if numel(unique(folders)) == 1
        % check if we took the whole folder
        if sum(tree{treeLevel}.clustering == unique(folders)) == numel(indsOfTile)
            isTileFitsByTree = true;
            return;
        end
    end
%    % merge
%    superFolders = tree{treeLevel}.super_folders(indsOfTile);
%    % if the last super folder is not all inside, we cannot merge 
%    if sum(tree{treeLevel+1}.clustering == max(superFolders)) ~=   sum(superFolders==max(superFolders))
%         return;
%    end

end
end

function currErr = evalTilingErr(data, tiling)
meanData = data;
for ti = 1:max(tiling.isbusy(:))
    [i_row, i_col] = find(tiling.isbusy==ti);
    meanData(i_row, i_col) = mean(mean(data(i_row, i_col)));
end
currErr = sum(sum((data - meanData).^2));
end
