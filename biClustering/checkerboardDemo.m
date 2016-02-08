
addpath(genpath('..\Questionnaire'));

n = 200;
p = 200;
mu = 2;
sigma = .01;
% Mu_rc = mu * rand(p, n);
[x,y]=meshgrid(linspace(-mu, mu, 10));
z = -x.*y + y;
Mu_rc = kron(z, ones(20));
X = randn(p, n)*sigma + Mu_rc;
if 0
row_perm = randperm(n);
col_perm = randperm(p);
else
   row_perm = 1:n;
col_perm = 1:p;
end
data = X(row_perm,:); 
data = data(:,col_perm);




%% COBRA
% MAXITERS = 100;
% tau = 1;
% U = data;
% Q = zeros(p, n);
% P = zeros(p, n);
% W = exp(-squareform(pdist(X.')));
% W = W/sum(W(:))/sqrt(n);
% W_tild = exp(-(squareform(pdist(X))));
% W_tild = W_tild/sum(W_tild(:))/sqrt(p);
% gam=1.5;
% for iters = 1:MAXITERS
% %     Y = prox_clust(U.' + P.', gam, W_tild).';
%     Y = prox_clust_l1(U.' + P.', gam);
%     P = U + P - Y;
% %     U = prox_clust(Y + Q, gam, W).';
%     U = prox_clust_l1(Y + Q, gam);
%     Q = Y + Q - U;
%     L(iters) = sum(sum((U - Y.').^2));
%     if L(iters) <= tau
%         break;
%     end
% end
%% Qu.

[ params ] = SetQuestParams( 0, 0, 0, 0 );
params.init_aff_row.on_rows = false;
params.init_aff_col = params.init_aff;
params.init_aff_col.metric = 'euc';
params.init_aff_row.metric = 'euc';
params.col_tree.eigs_num = 12;
params.row_tree.eigs_num = 12;
params.n_iters = 2;
[ row_tree, col_tree ] = RunQuestionnaire( params, data );
row_aff = CalcEmdAff(data.', col_tree, params.row_emd);
col_aff = CalcEmdAff(data, row_tree, params.col_emd);

[ err_rate ] = OrganizeData( X, data, row_aff, col_aff, row_perm, col_perm, 2, 2  );


