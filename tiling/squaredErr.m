function representation_error = squaredErr(data, est)

diff2 = (data-est).^2;
summing = sum(diff2(:));
representation_error = (summing);