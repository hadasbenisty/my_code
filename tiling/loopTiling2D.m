function [minErr, solutionTiling] = loopTiling2D(data, row_tree, col_tree, vol)


minErr = Inf;solutionTiling=[];
MAX_S = numel(data);
MAX_ITERS = 1e6;
stind = 1;
stck(stind).tiling.isbusy = zeros(size(data, 1), size(data, 2), 1);
stck(stind).tiling.isLeader = zeros(size(data, 1), size(data, 2), 1);
stck(stind).ind2data = 1;
tilesList = getTiles(vol);
stck(stind).ti = 1;

for iter = 1:MAX_ITERS
    
    if stck(stind).ind2data > MAX_S ||  sum(stck(stind).tiling.isbusy(:) == 0) == 0
        stck(stind).currErr = evalTilingErr(data, stck(stind).tiling);
        if stck(stind).currErr < minErr
            minErr =  stck(stind).currErr;
            disp(['Volume = ' num2str(vol) ' err = ' num2str(minErr)]);
            %         disp(tiling.isbusy);
            disp('---------------------');
            solutionTiling = stck(stind).tiling;
        end
        stind = stind - 1;
        stck(stind).tiling = unmarkTile(stck(stind).tiling, tilesList(stck(stind).ti), stck(stind).ind2data);
        stck(stind).ti =  stck(stind).ti + 1;
        
    end
    [row_i,col_i] = ind2sub(size(data),stck(stind).ind2data);
    if stck(stind).tiling.isbusy(row_i, col_i) % check if it's busy and shift to the next one
        stck(stind).ind2data = stck(stind).ind2data + 1;
        stck(stind).ti = 1;
        continue;
    end
    
    
    while stck(stind).ti <= length(tilesList)
        % cheking if this tile fits
        isTileFits = checkTile(data, row_tree, col_tree, stck(stind).ind2data, tilesList(stck(stind).ti));
        if isTileFits
            % marking this tile
            stck(stind).tiling = markTile(stck(stind).tiling, tilesList(stck(stind).ti), stck(stind).ind2data);
            currErr = evalTilingErr(data, stck(stind).tiling);
            if currErr < minErr
                stck(stind+1) = stck(stind);
                stck(stind+1).ind2data = stck(stind+1).ind2data + 1;
                stck(stind+1).ti = 1;
                stind = stind + 1;
                break;
            else
                stck(stind).tiling = unmarkTile(stck(stind).tiling, tilesList(stck(stind).ti), stck(stind).ind2data);
            end
        end
        stck(stind).ti = stck(stind).ti + 1;
    end
    if stck(stind).ti > length(tilesList)
        if stind == 1
            return;
        end
        % remove this tile and continue to examine the other tiles in the
        % list
        stind = stind - 1;
        stck(stind).tiling = unmarkTile(stck(stind).tiling, tilesList(stck(stind).ti), stck(stind).ind2data);
        stck(stind).ti =  stck(stind).ti + 1;
    end
end






