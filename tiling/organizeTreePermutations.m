function orederedTree = organizeTreePermutations(origData, X, Tree, treeLevel, err_rate)

if isbad(err_rate)
    err_rate = CalcErrRate( origData, X );
end
% we have reached the highest level of the tree - are finished
if Tree{treeLevel}.folder_count == 1
    orederedTree = Tree;
    return;
end

% take all super folders
orederedTree = Tree;
sups = unique(Tree{treeLevel}.super_folders);
for supi = 1:length(sups)
    % inds to the clusters we are going to switch
   inds4sups = find( Tree{treeLevel}.super_folders == sups(supi));
   perms4sups = perms(1:length(inds4sups));
   for permi = 1:size(perms4sups, 1)
     curr_orederedTree = permTree(orederedTree, treeLevel, perms4sups(permi, :));
     [~, currclustering] = sort(curr_orederedTree{2}.clustering, 'ascend');
     reorderedData = X(currclustering, :);
     curr_err_rate  = CalcErrRate( origData, reorderedData );
     % saving this switch
     if curr_err_rate < err_rate
         err_rate = curr_err_rate;
         orederedTree = curr_orederedTree;
     end  
   end
end
