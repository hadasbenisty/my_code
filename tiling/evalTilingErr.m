function currErr = evalTilingErr(data, tiling)
meanData = data;
for ti = 1:max(tiling.isbusy(:))
    [i_row, i_col] = find(tiling.isbusy==ti);
    meanData(i_row, i_col) = mean(mean(data(i_row, i_col)));
end
currErr = sum(sum((data - meanData).^2));
end