function plotSlices(X, numOfSlices, colttl, rowttl, trialttl)

if nargin == 1
    numOfSlices = 4;
end
LenDims = length(size(X));
if LenDims == 3
      [n_cols, n_rows, nTrials] = size(X);
   
        pltslices(permute(X, [1 2 3 ]), numOfSlices, colttl, rowttl, trialttl);
        pltslices(permute(X, [2 3 1 ]), numOfSlices, rowttl, trialttl, colttl);
        pltslices(permute(X, [3 1 2 ]), numOfSlices, trialttl, colttl, rowttl);
 
else
    [dims, n_cols, n_rows, nTrials] = size(X);
    for d = 1:dims
        pltslices(permute(X(d, :, :, :), [2 3 4 1]), numOfSlices, colttl, rowttl, trialttl);
        
    end
    
    for d = 1:dims
        pltslices(permute(X(d, :, :, :), [3 4 2 1]), numOfSlices, rowttl, trialttl, colttl);
    end
    
    for d = 1:dims
        pltslices(permute(X(d, :, :, :), [4 2 3 1]), numOfSlices, trialttl, colttl, rowttl);
    end
end

end
