function [ind_H, ind_V, ind_D] = get_subtrees(m)

% number of levels
L = log2(m^2-1)/2;

n_subtree = (m^2-1)/3;
ind_H = zeros(n_subtree,2);
ind_V = zeros(n_subtree,2);
ind_D = zeros(n_subtree,2);

ind_H(1,:) = [1, 2];
ind_V(1,:) = [2, 1];
ind_D(1,:) = [2, 2];

cont = 1;
for ell = 1:L
   % Number of nodes
   N = 4^(ell-1);
   % Nodes indices
   ind_nodes = (1:N)+(4^(ell-1)-1)/3;
   for n = 1:N      
      temp = coeff_H_children(ind_H(ind_nodes(n),:));
      n_children = size(temp,1);
      ind_H(cont+1:cont+n_children,:) = temp;
      temp = coeff_H_children(ind_V(ind_nodes(n),:));
      ind_V(cont+1:cont+n_children,:) = temp;
      temp = coeff_H_children(ind_D(ind_nodes(n),:));
      ind_D(cont+1:cont+n_children,:) = temp;
      cont = cont + n_children;
   end
end

function ind = coeff_H_children(ind)

x = ind(1);
y = ind(2);

ax = 2*x-1;
ay = 2*y-1;
ind(1,:) = [ax, ay];
ind(2,:) = [ax, ay+1];
ind(3,:) = [ax+1, ay];
ind(4,:) = [ax+1, ay+1];