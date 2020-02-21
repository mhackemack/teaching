function out = binomial(i,n)
f = factorial([n,i,n-i]);
out = f(1) / f(2) / f(3);