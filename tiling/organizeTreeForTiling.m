function orderedtree = organizeTreeForTiling(Tree, ordering)

for treeLevel = 1:length(Tree)
    orderedtree{treeLevel} = Tree{treeLevel};
    orderedtree{treeLevel}.clustering = Tree{treeLevel}.clustering( ordering);
end