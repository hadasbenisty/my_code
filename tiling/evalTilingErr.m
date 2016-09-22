function [currErr, meanData] = evalTilingErr(data, tiling, params)
meanData = zeros(size(data));
w_mat = ones(size(data));
[row_sz, col_sz] = size(data);
tiles = unique(tiling.isbusy(:));
tiles = setdiff(tiles, 0);
for ti = 1:length(tiles)
    [i_row, i_col] = find(tiling.isbusy==tiles(ti));
    w_row = ((max(i_row)-min(i_row) + 1)/row_sz)^(params.beta_row );
    w_col = ((max(i_col)-min(i_col) + 1)/col_sz)^(params.beta_col);
    meanData(min(i_row):max(i_row), min(i_col):max(i_col)) = mean(mean(data(min(i_row):max(i_row), min(i_col):max(i_col))));
    w_mat(min(i_row):max(i_row), min(i_col):max(i_col)) = w_row*w_col;
end
% currErr = sum(sum(((data(nnzero_inds) - meanData(nnzero_inds)).^2).*w_mat(nnzero_inds)));
currErr = sum(sum((abs(data(meanData~=0) - meanData(meanData~=0))).*w_mat(meanData~=0)));

if isbad(currErr)
    error('Something is wrong with the error calculation');
end
% currErr1 = sum(sum((abs(data - meanData))./denum_mat));

% currErr = sum(sum(((data - meanData).^2)));

end