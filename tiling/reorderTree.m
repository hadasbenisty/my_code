function orderedtree = reorderTree(tree, order)

orderedtree = tree;
for tree_level = length(orderedtree)-1:-1:1
    curr_super_folders = orderedtree{tree_level}.super_folders;
    for sfind = 1:orderedtree{tree_level+1}.folder_count
        currinds = find(curr_super_folders==sfind);
        [~, ordered_clusters] = sort(order(currinds), 'descend');
        orderedtree{tree_level}.super_folders(currinds) = orderedtree{tree_level}.super_folders(currinds(ordered_clusters));
        
        orderedtree{tree_level}.folder_sizes(currinds) = orderedtree{tree_level}.folder_sizes(currinds(ordered_clusters));
    end
    orderedtree{tree_level}.folder_count = orderedtree{tree_level}.folder_count;
    orderedtree{tree_level}.clustering = orderedtree{tree_level}.clustering;
    
end
orderedtree{ length(tree)} = tree{end};
tree1{2} = tree{end};
tree1(1) = orderedtree(3);

subplot(2,1,1);treeplot(nodes(tree),'.')
subplot(2,1,2);treeplot(nodes(tree1),'.')