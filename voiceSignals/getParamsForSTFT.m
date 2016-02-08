function params = getParamsForSTFT(eigsnum_col,eigsnum_row, eigsnum_trials)

params.verbose = 2;
params.init_aff_col.on_rows = false;
params.init_aff_col.metric = 'cosine_similarity';
% params.init_aff_col.metric = 'Euc';
params.init_aff_col.eps = 2;
params.init_aff_col.knn = 10;
params.init_aff_col.thresh = 0.0;

params.init_aff_row.on_rows = false;
% params.init_aff_row.metric = 'cosine_similarity';
params.init_aff_row.metric = 'Euc';
params.init_aff_row.eps = 1;
params.init_aff_row.knn = 50;
params.init_aff_row.thresh = 0;

params.init_aff_trials.on_rows = false;
params.init_aff_trials.metric = 'cosine_similarity';
% params.init_aff_trials.metric = 'Euc';
params.init_aff_trials.eps = 1;
params.init_aff_trials.knn = 50;
params.init_aff_trials.thresh = 0;
params.row_tree.verbose = 0;

params.row_tree.eigs_diff_th = 0.05;
params.row_tree.constant=2;
params.row_tree.min_joins_percentage = 0.4;
params.row_tree.eigs_num = eigsnum_row;
params.col_tree.constant=1;
params.col_tree.eigs_num = eigsnum_col;
params.col_tree.min_joins_percentage = 0.5;
params.col_tree.eigs_diff_th = 0.05;
params.col_tree.verbose = 0;

params.trials_tree.constant=2;
params.trials_tree.eigs_num = eigsnum_trials;
params.trials_tree.min_joins_percentage = 0.5;
params.trials_tree.eigs_diff_th = 0.05;
params.trials_tree.verbose = 0;

params.n_iters = 3;
params.col_emd.beta=0;
params.col_emd.alpha = 0;
params.col_emd.eps = 2;
params.trials_emd.beta=0;
params.trials_emd.alpha = 0;
params.trials_emd.eps = 10;

params.row_emd.beta=0;
params.row_emd.alpha = 1;
params.row_emd.eps = 2;
params.data.to_normalize = true;
params.data.over_rows = true;
params.data.normalization_type = 'by_std';