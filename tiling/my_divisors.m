function D = my_divisors(n)


K=1:n;
D = K(rem(n,K)==0);