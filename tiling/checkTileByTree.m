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
