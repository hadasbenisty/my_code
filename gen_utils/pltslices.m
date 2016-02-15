function pltslices(X, numOfSlices, colttl, rowttl, trialttl)

n_cols = size(X, 1);

inds = round(linspace(1, n_cols, numOfSlices));

   slc = X(inds, :, :);  
   R = ceil(sqrt(numOfSlices));
   
figure;
for n=1:numOfSlices
   subplot(R,  R, n)
   imagesc(shiftdim(slc(n, :, : )));colorbar;title(['Slice no ' num2str(n) ' Of ' colttl]);
   ylabel(rowttl);
   xlabel(trialttl);
end
end