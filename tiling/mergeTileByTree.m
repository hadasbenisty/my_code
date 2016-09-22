function new_r = mergeTileByTree(tree, inds_r)

for l=2:length(  tree)
    cl = tree{l}.clustering(inds_r) ;
    if numel(unique(cl)) ==1
        break;
    end
end
for l1 = l + 1: length(tree)
    new_r =find(tree{l1}.clustering == cl(1));
    if numel(new_r) > inds_r
        break;
    end
end
end