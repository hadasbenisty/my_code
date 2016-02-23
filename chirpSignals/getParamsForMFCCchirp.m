function params = getParamsForMFCCchirp(eigsnum_col,eigsnum_row)

    params.init_aff.on_rows = false;
    params.init_aff.metric = 'cosine_similarity';
    params.init_aff.metric = 'Euc';
    params.init_aff.eps = 4;
    params.init_aff.knn = 10;
    params.init_aff.thresh = 0;
   
    
    %mbm,b,mnmk
    params.row_tree.eigs_num = 75;
    params.row_tree.constant=2;
    params.row_tree.min_joins_percentage = 0.4;
    params.row_tree.eigs_num = eigsnum_row;
    params.col_tree.constant=1;
    params.col_tree.eigs_num = eigsnum_col;
    params.col_tree.min_joins_percentage = 0.5;
    params.n_iters = 4;
    params.col_emd.beta=0;
    params.col_emd.alpha = 1;
    params.col_emd.eps = 2;
    params.row_emd.beta=0;
    params.row_emd.alpha = 1;
    params.row_emd.eps = 2;
    params.data.to_normalize = true;
    params.data.over_rows = true;
    params.data.normalization_type = 'by_std';