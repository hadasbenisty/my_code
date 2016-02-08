function om = omega_fun(W, U)

om = 0;
for j = 1: size(U, 2)
    for i = 1: j-1
    om = om + W(i, j) * sqrt(sum((U(:, i) - U(:, j)).^2));
    end
end