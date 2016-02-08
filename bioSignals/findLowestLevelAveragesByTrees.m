function U = findLowestLevelAveragesByTrees(col_tree, row_tree, X)

U = zeros(size(X));

for nclusters_cols = 1:col_tree{1}.folder_count
    for nclusters_rows = 1:row_tree{1}.folder_count
        U(row_tree{1}.clustering == nclusters_rows, col_tree{1}.clustering == nclusters_cols) = X(row_tree{1}.clustering == nclusters_rows, col_tree{1}.clustering == nclusters_cols);
    end
end
    
