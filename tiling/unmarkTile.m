function tiling = unmarkTile(tiling, tile, ind2data)
% disp(['Unmarking tile ' num2str(tile.row) ' by ' num2str(tile.col) ' in ' num2str(ind2data)]);

[row_i,col_i] = ind2sub(size(tiling.isbusy),ind2data);
tiling.isLeader(row_i,col_i) = false;
tiling.isbusy(row_i:row_i+tile.row-1, col_i:col_i+tile.col-1) = 0;
imagesc(tiling.isbusy);title(['Tiling Vol = ' num2str(tile.row*tile.col)]);
drawnow;

end
