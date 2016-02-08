function V = prox_clust_l1(Z, sig)

V = max(1-sig./(abs(Z) + eps), 0).*Z;