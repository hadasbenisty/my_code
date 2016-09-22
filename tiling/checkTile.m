function isTileFits = checkTile(data, row_tree, col_tree, ind2data, tile)

% first assuming false
isTileFits = false;

[size_row, size_col, ~] = size(data);
[row_i,col_i] = ind2sub([size_row, size_col],ind2data);

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
