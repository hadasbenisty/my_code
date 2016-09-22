function [meanMat, allMat, meanMatAlltrials] = getClusteringByTreeAndPlot(Trees, data, NeuronsLabels, newNeuronsLabels, params, ttl)

[meanMat, allMat, meanMatAlltrials] = getCentroidsByTree(Trees, data, NeuronsLabels, newNeuronsLabels);
[aff_mat  ] = CalcInitAff( meanMat.', params );
[vecs, vals] = CalcEigs(threshold(aff_mat, 0.0), 2);%
[ row_order ] = OrganizeDiffusion3DbyOneDim( meanMat, vecs*vals );

plotByClustering(meanMat(row_order, :),  ttl);
plotByClustering(allMat(row_order),  ttl);