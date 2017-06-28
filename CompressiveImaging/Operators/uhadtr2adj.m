function y = uhadtr2adj(x, ind, N)

y = zeros([N N]);
y(ind) = x;
y = uhadtr2(y);
