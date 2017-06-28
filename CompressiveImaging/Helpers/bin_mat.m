function B = bin_mat(d,m,n)

nonZeros = d;

rows = zeros(n*nonZeros,1);
cols = zeros(n*nonZeros,1);
vals = true(n*nonZeros,1);

tic
for j=1:n
  cols(1+(j-1)*nonZeros:j*nonZeros) = j*ones(nonZeros,1);
  rows(1+(j-1)*nonZeros:j*nonZeros) = randperm(m,nonZeros);
  if mod(j,1e5) == 0
     fprintf('Expander completion %.2f %%\n',j/n*100);
  end
end
toc

B = sparse(rows,cols,vals,m,n);
