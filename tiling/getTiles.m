function tilesList = getTiles(vol)


D = my_divisors(vol);

for k = 1:length(D)
    tilesList(k).row = D(k);
    tilesList(k).col = vol/tilesList(k).row;
end