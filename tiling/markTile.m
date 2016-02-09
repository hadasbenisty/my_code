function tiling = markTile(tiling, tile, ind2data)
       
%          /* put tile on ind2data */
%          mark ind2data as leader of tile;
disp(['Marking tile ' num2str(tile.row) ' by ' num2str(tile.col) ' in ' num2str(ind2data)]);
[row_i,col_i] = ind2sub(size(tiling.isbusy),ind2data);
tiling.isLeader(row_i,col_i) = true;
tiling.isbusy(row_i:row_i+tile.row-1, col_i:col_i+tile.col-1) = max(tiling.isbusy(:))+1;
imagesc(tiling.isbusy);
title(['Tiling Vol = ' num2str(tile.row*tile.col)]);
drawnow;

%          mark all neighbors M[sx .. sx+tx-1][sy .. sy+ty-1] = busy;
end
