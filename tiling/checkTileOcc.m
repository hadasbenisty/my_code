function isTileNotOcc = checkTileOcc(tiling, ind2data, tile)


[row_i,col_i] = ind2sub(size(tiling.isbusy),ind2data);
if any(tiling.isbusy(row_i:row_i+tile.row-1, col_i:col_i+tile.col-1)~= 0)
    isTileNotOcc = false;
else
    isTileNotOcc = true;
end