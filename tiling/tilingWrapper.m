function [meanTiled, solutionTiling, n1_roder, n2_roder, data4tilingLog] = tilingWrapper(representation_error_fun, verbose, Tree1, Tree2, eps, x)
p1.representation_error_fun = representation_error_fun;
p1.verbose = verbose;
[~, n1_roder] = sort(Tree1{2}.clustering);
n2_roder = getTimeOrderFromTree(Tree2{2}.clustering);
% [~, n3_roder] = sort(Trees{runningOrder==3}{2}.clustering);

n1_orderedtree = organizeTreeForTiling(Tree1, n1_roder);
n2_orderedtree = organizeTreeForTiling(Tree2, n2_roder);


p1.eps=eps;
orderedData = x(n1_roder, :, :);%
orderedData = orderedData(:, n2_roder, :);
data4tiling = mean(orderedData,3);
data4tilingLog = log(threshold(data4tiling+2*abs(min(data4tiling(:))),0)+eps);

if p1.verbose>1
    figure;
    subplot(2,2,3);plotTreeWithColors(Trees{runningOrder==1}, numLabels(n1_roder));title('wires');view(-90,90);
    subplot(2,2,2);plotTreeWithColors(Trees{runningOrder==2}, 1:62);title('Time');
    subplot(2,2,4);
end
[~, solutionTiling] = reversedIterativeTilingGenericDimsleagal(data4tilingLog, {n1_orderedtree, n2_orderedtree}, p1);
[meanTiled] = getTiledData(data4tiling, solutionTiling);
