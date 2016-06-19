function [minErr, tiling, solutionTiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data, tiling, solutionTiling, minErr)

params.beta_row=0;
params.beta_col=0;
MAX_S = numel(data);
if ind2data > MAX_S ||  sum(tiling.isbusy(:) == 0) == 0
    currErr = evalTilingErr(data, tiling, params);
    if currErr < minErr
        minErr =  currErr;
%         disp(['Volume = ' num2str(vol) ' err = ' num2str(minErr)]);
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
% tilesList(1).row=10;
% tilesList(1).col=70;
% tilesList(2).row=20;
% tilesList(2).col=30;
% tilesList(3).row=60;
% tilesList(3).col=10;
% 
% tilesList(4).row=40;
% tilesList(4).col=40;
% 
% tilesList(5).row=40;
% tilesList(5).col=20;
% 
% tilesList(6).row=20;
% tilesList(6).col=20;
% 
% tilesList(7).col=20;
% tilesList(7).col=20;
% tilesList=tilesList(1:7);
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
        tiling = markTile(tiling, tilesList(ti), ind2data, 2);
        % see if we don't already get an error higher than the prev. so
        % there is a point to continue
        currErr = evalTilingErr(data, tiling, params);
        if currErr < minErr
            [minErr, tiling, solutionTiling] = recursiveTiling2D(data, row_tree, col_tree, vol, ind2data+1, tiling, solutionTiling, minErr);
        end
        % remove this tile and continue to examine the other tiles in the
        % list
        tiling = unmarkTile(tiling, tilesList(ti), ind2data, 2);
    end
end
end






