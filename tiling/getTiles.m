function tilesList = getTiles(vol)

l = 1;
for vi = 1:length(vol)
    D = my_divisors(vol(vi));

for k = 1:length(D)
    tilesList(l).row = D(k);
    tilesList(l).col = vol(vi)/tilesList(l).row;
    l = l + 1;
end
end