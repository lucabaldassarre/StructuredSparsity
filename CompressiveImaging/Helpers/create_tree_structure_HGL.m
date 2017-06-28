function tree = create_tree_structure_HGL(m, D)

N = m^2;
n_subtree = (m^2-1)/3;
depth = log(n_subtree*(D-1)+1)/log(D) - 1;

tree.eta_g = ones(N,1);
tree.N_own_variables = int32(ones(N,1));
tic
[ind_H, ind_V, ind_D] = get_subtrees(m);
toc
tic
ind_H = sub2ind([m,m],ind_H(:,1), ind_H(:,2));
ind_V = sub2ind([m,m],ind_V(:,1), ind_V(:,2));
ind_D = sub2ind([m,m],ind_D(:,1), ind_D(:,2));
toc
tree.own_variables = int32([1; ind_H; ind_V; ind_D])-1;

tic
i = [2; 2+n_subtree; 2+2*n_subtree];
j = [1; 1; 1];

i = [i; (3:n_subtree+1)'];
jt = 2:(n_subtree+1-4^depth);
jt = repmat(jt,4,1);
jt = jt(:);
j = [j; jt];

i = [i; (n_subtree+3:2*n_subtree+1)'];
jt = n_subtree+2:2*n_subtree+1-4^depth;
jt = repmat(jt,4,1);
jt = jt(:);
j = [j; jt];

i = [i; (2*n_subtree+3:3*n_subtree+1)'];
jt = 2*n_subtree+2:3*n_subtree+1-4^depth;
jt = repmat(jt,4,1);
jt = jt(:);
j = [j; jt];

s = ones(numel(i),1);

tree.groups = sparse(i,j,s,n_subtree*3+1,n_subtree*3+1);
toc
