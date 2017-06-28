function y = dct2adj(x, ind, N)

y = zeros([N N]);
y(ind) = x;
y = idct2(y);
