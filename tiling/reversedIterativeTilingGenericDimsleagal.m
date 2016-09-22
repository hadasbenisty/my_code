%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A tiling algorithm based on trees J.B, with top down and genric dims.
% Searches for the "best" tiling w.r.t. a cost:
%              cost = representation_error + eps*number_of_tiles
% inputs:
% data   - samples, of any dimmension 
% Trees  - a cell array such that length(Trees) == length(size(data))
% params - parameters struct with the following fields:
%            * max_tiles_num            - maximal amount of tiles (default - Inf)
%            * eps                      - weight for the amount of tiles in the cost
%            * representation_error_fun - a handle for a function
%            evaluating the representation error (||_1, or ||_2 for
%            example)
% 
% written by: Hadas Benisty, hadas.benisyu@gmail.com 
% v. 1.0 : 19/6/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [currCost, tilingSol, allsols] = reversedIterativeTilingGenericDimsleagal(data, Trees, params)

if ~isfield(params, 'max_tiles_num')
    params.max_tiles_num = Inf;
end
isfinished = false;
tilingSol.isbusy = ones(size(data));
allsols(1).isbusy = tilingSol.isbusy;
currCost = evalCost(data, tilingSol.isbusy, params.representation_error_fun, params.eps);
MAX_ITERS = 100;

for iters = 1:MAX_ITERS
    disp(iters);
    tiles = unique(tilingSol.isbusy(:));
    
    splitting = zeros(length(tiles), 1);
    for ti = 1:length(tiles)
        sizeTiles(ti) = sum(tilingSol.isbusy(:) == tiles(ti));
    end
    [~,indsSort] = sort(sizeTiles, 'descend');
%     tiles=tiles(indsSort);
    for ti = 1:length(tiles)
        %         [inds_r, inds_c] = find(tilingSol.isbusy == tiles(ti));
        if isempty(find(tilingSol.isbusy==tiles(ti)))
            warning('1');
            continue;
        end
        splits = getNextSplitFromTree(Trees, tilingSol.isbusy, tiles(ti));
        cost = zeros(length(splits), 1);
        for n = 1:length(splits)
            if length(unique(splits{n})) == 1
                cost(n) = inf;
            else
                possibleTilingSol(n).isbusy = splitTile(tilingSol.isbusy, tiles(ti), splits, n);
                cost(n) = evalCost(data, possibleTilingSol(n).isbusy, params.representation_error_fun, params.eps);
            end
        end
        [~, mini] = min([currCost; cost]);
        if mini == 1
            % do nothing
        else
            currCost = cost(mini-1);
            tilingSol = possibleTilingSol(mini-1);
            allsols(end+1).isbusy = tilingSol.isbusy;
            splitting(ti) = 1;
            disp(currCost);
            if params.verbose > 1
                meanTiled = getTiledData(data, tilingSol);
                imagesc(meanTiled(:,:,1));
                drawnow;
            end
            if length(unique(tilingSol.isbusy(:))) >= params.max_tiles_num
                isfinished = true;
                break;
            end
        end
    end
    if all(splitting == 0) || isfinished
        break;
    end
end
end
function splits = getNextSplitFromTree(Trees, isbusy, tile)

splits = cell(length(Trees), 1);
inds = find(isbusy == tile); %#ok<NASGU>
cmdStr = '[';
for n = 1:length(Trees)
    cmdStr = [cmdStr, 'subinds' num2str(n) ',']; %#ok<AGROW>
end
cmdStr = [cmdStr(1:end-1), '] = ind2sub(size(isbusy), inds);'];
eval(cmdStr);
finaltreeLevel = zeros(length(Trees), 1);
for n = 1:length(Trees)
    currsubinds = eval(['subinds' num2str(n)]);
    len(n) = length(unique(currsubinds));
    
    for tree_level = 1:length(Trees{n})
        if length(unique(Trees{n}{tree_level}.clustering(currsubinds(1):currsubinds(1) + len(n) -1))) == 1
            finaltreeLevel(n) = tree_level;
            break;
        end
    end
end
for n = 1:length(Trees)
    currsubinds = eval(['subinds' num2str(n)]);
    if finaltreeLevel(n) == 1
        splits{n} = 1;
    else
        splits{n} = Trees{n}{finaltreeLevel(n)-1}.clustering(currsubinds(1):currsubinds(1) + len(n) -1);
    end
end


end

function newIsbusy = splitTile(isbusy, tile, splits, n2split)

newsplits = cell(size(splits));
for n = 1:length(splits)
    newsplits{n} = ones(size(splits{n}));
end
newsplits{n2split} = splits{n2split};
newIsbusy = isbusy;
inds = find(isbusy == tile); %#ok<NASGU>
cmdStr = '[';
for n = 1:length(newsplits)
    cmdStr = [cmdStr, 'subinds' num2str(n) ',']; %#ok<AGROW>
end
cmdStr = [cmdStr(1:end-1), '] = ind2sub(size(isbusy), inds);'];
eval(cmdStr);

for n = 1:length(newsplits)
    currsubinds = eval(['subinds' num2str(n)]);
    tile_inds{n} = min(currsubinds):max(currsubinds);
    clusters{n} = unique(newsplits{n});
end
cmdStr = ['for ci'];
for n = 1:length(newsplits)
    cmdStr = [cmdStr num2str(n) ' = 1:length(clusters{' num2str(n) '}), for ci'];
end
cmdStr = cmdStr(1:end-6);
cmdStr = [cmdStr 'newIsbusy(tile_inds{'];
for n = 1:length(newsplits)
    cmdStr=[cmdStr num2str(n) '}(find(newsplits{' num2str(n) '} == clusters{' num2str(n) '}(ci' num2str(n) '))), tile_inds{'];
end
cmdStr=[cmdStr(1:end-12) ')= max(newIsbusy(:)) + 1;'];
for n = 1:length(newsplits)
    cmdStr = [cmdStr ' end; '];
end
eval(cmdStr);
end
