function [meanData] = getTiledData(data, tiling)
meanData = zeros(size(data));
tiles = unique(tiling.isbusy(:));
tiles = setdiff(tiles, 0);
for ti = 1:length(tiles)
    inds = find(tiling.isbusy==tiles(ti));
    cmdStr = '[';
    for n = 1:length(size(data))
        cmdStr = [cmdStr, 'subinds' num2str(n) ',']; %#ok<AGROW>
    end
    cmdStr = [cmdStr(1:end-1), '] = ind2sub(size(tiling.isbusy), inds);'];
    eval(cmdStr);
    
    
    cmdStr = 'meanData(min(subinds';
    for n = 1:length(size(data))
        cmdStr=[cmdStr num2str(n) '):max(subinds' num2str(n) '), min(subinds'];
    end
    cmdStr=[cmdStr(1:end-13) ') = '];
    for n = 1:length(size(data))
        cmdStr=[cmdStr 'mean('];
    end
    cmdStr=[cmdStr 'data(min(subinds'];
    for n = 1:length(size(data))
        cmdStr=[cmdStr num2str(n) '):max(subinds' num2str(n) '), min(subinds'];
    end
    cmdStr=[cmdStr(1:end-13) ')'];
for n = 1:length(size(data))
cmdStr=[cmdStr ')'];
end
cmdStr=[cmdStr ';'];
    eval(cmdStr);
end

end