function params = getParamsForBioData(eigsnum_col,eigsnum_row, eigsnum_trials)
params.init_aff_col.on_rows = false;
params.init_aff_col.metric = 'cosine_similarity';
% params.init_aff_col.metric = 'Euc';
params.init_aff_col.eps = 2;
params.init_aff_col.knn = 10;
params.init_aff_col.thresh = 0.2;

params.init_aff_row.on_rows = false;
params.init_aff_row.metric = 'cosine_similarity';
params.init_aff_row.eps = 2;
params.init_aff_row.knn = 10;
params.init_aff_row.thresh = 0;
params.verbose = 2;
params.col_tree.verbose = 2;
params.row_tree.verbose = 2;
params.trials_tree.verbose = 2;
params.init_aff_trials.on_rows = false;
params.init_aff_trials.metric = 'cosine_similarity';
params.init_aff_trials.eps = 2;
params.init_aff_trials.knn = 10;
params.init_aff_trials.thresh = 0;

params.row_tree.constant=2;
params.row_tree.min_joins_percentage = 0.2;
params.row_tree.eigs_num = eigsnum_row;
params.col_tree.constant=1;
params.col_tree.eigs_num = eigsnum_col;
params.col_tree.min_joins_percentage = 0.2;

params.trials_tree.constant=2;
params.trials_tree.eigs_num = eigsnum_trials;
params.trials_tree.min_joins_percentage = 0.2;

params.n_iters = 2;
params.col_emd.beta=0;
params.col_emd.alpha = 0;
params.col_emd.eps = 2;
params.trials_emd.beta=0;
params.trials_emd.alpha = 0;
params.trials_emd.eps = 1;

params.row_emd.beta=0;
params.row_emd.alpha = 0;
params.row_emd.eps = 2;
params.data.to_normalize = false;
params.data.over_rows = true;
params.data.normalization_type = 'by_std';