function params = getParamsForSTFTchirp(eigsnum_col,eigsnum_row, eigsnum_trials)

% params.verbose = 0; % no disp and no plot
% params.verbose = 1; % disp and no plot
params.verbose = 2; % disp and plot
params.init_aff_col.on_rows = true;
params.init_aff_col.metric = 'cosine_similarity';
params.init_aff_col.metric = 'Euc';
params.init_aff_col.eps = 5;
params.init_aff_col.knn = 50;
params.init_aff_col.thresh = 0.3;%-inf;%0.6;%0.3;

params.init_aff_row.on_rows = false;
params.init_aff_row.metric = 'cosine_similarity';
% params.init_aff_row.metric = 'Euc';
params.init_aff_row.eps = 1;
params.init_aff_row.knn = 50;
params.init_aff_row.thresh = 0.3;%-0.2;%-inf;% -0.2;

params.init_aff_trials.on_rows = false;
params.init_aff_trials.metric = 'cosine_similarity';
% params.init_aff_trials.metric = 'Euc';
params.init_aff_trials.eps = 1;
params.init_aff_trials.knn = 50;
params.init_aff_trials.thresh = 0.3;%-0.2;%-inf;% -0.2;

params.row_tree.eigs_diff_th = 0.05;

params.row_tree.verbose = 0;

params.row_tree.constant=0.5;%2;
params.row_tree.min_joins_percentage = 0.5;%0.1
params.row_tree.eigs_num = eigsnum_row;
params.col_tree.verbose = 0;

params.col_tree.constant=0.5;%1;
params.col_tree.eigs_num = eigsnum_col;
params.col_tree.min_joins_percentage = 0.5;%0.1;
params.col_tree.eigs_diff_th = 0.05;
params.trials_tree.eigs_diff_th = 0.05;
params.trials_tree.verbose = 0;

params.trials_tree.constant=2;
params.trials_tree.eigs_num = eigsnum_trials;
params.trials_tree.min_joins_percentage = 0.1;%.50

params.n_iters = 2;
params.col_emd.beta=0;
params.col_emd.alpha = 1;%0 works well for 2D
params.col_emd.eps = .5;
params.trials_emd.beta=0;
params.trials_emd.alpha = 1;%0 works well for 2D
params.trials_emd.eps = 1;

params.row_emd.beta=0;
params.row_emd.alpha =0; -1;
params.row_emd.eps = .5;
params.data.to_normalize = true;
params.data.over_rows = false;
params.data.normalization_type = 'by_std';