function     loc = findStrLocInCellArray(cellA, strList)

loc = zeros(length(strList), 1);
for k = 1:length(strList)
    currLoc = find(strcmp(strList{k}, cellA));
    if ~isempty(currLoc)
        loc(k) = currLoc;
    end
end