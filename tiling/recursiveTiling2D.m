function [minErr, tiling, solutionTiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data, tiling, solutionTiling, minErr)


MAX_S = numel(data);
if ind2data > MAX_S ||  sum(tiling.isbusy(:) == 0) == 0
    currErr = evalTilingErr(data, tiling);
    if currErr < minErr
        minErr =  currErr;
        disp(['Volume = ' num2str(vol) ' err = ' num2str(minErr)]);
%         disp(tiling.isbusy);    
        disp('---------------------');
        solutionTiling = tiling;
        return;
    end
end
[row_i,col_i] = ind2sub(size(data),ind2data);
if tiling.isbusy(row_i, col_i) % check if it's busy and shift to the next one
    [minErr, tiling, solutionTiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data+1, tiling, solutionTiling, minErr);
    return
end

tilesList = getTiles(vol);
for ti = 1:length(tilesList)
    % cheking if this tile fits
    isTileFits = checkTile(data, row_tree, col_tree, ind2data, tilesList(ti));
    if ~isTileFits
        continue;
    end
    % checking if it is not occupied
    isTileNotOcc = checkTileOcc(tiling, ind2data, tilesList(ti));
    if isTileNotOcc
        % marking this tile 
        tiling = markTile(tiling, tilesList(ti), ind2data);
        % see if we don't already get an error higher than the prev. so
        % there is a point to continue
        currErr = evalTilingErr(data, tiling);
        if currErr < minErr
            [minErr, tiling, solutionTiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data+1, tiling, solutionTiling, minErr);
        end
        % remove this tile and continue to examine the other tiles in the
        % list
        tiling = unmarkTile(tiling, tilesList(ti), ind2data);
    end
end
end






