function [xpermed, n1_perm, n2_perm, n3_perm] = permdata(x, ispermute)

if ispermute
    n1_perm = randperm(size(x,1));
    n2_perm = randperm(size(x,2));
    n3_perm = randperm(size(x,3));
else
    n1_perm = 1:size(x,1);
    n2_perm = 1:size(x,2);
    n3_perm = 1:size(x,3);
end

xpermed = x(n1_perm,:,:);
xpermed = xpermed(:,n2_perm,:);
xpermed = xpermed(:,:,n3_perm);
