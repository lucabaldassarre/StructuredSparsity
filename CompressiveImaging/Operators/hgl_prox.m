function out = hgl_prox(x, t, tree, param)

param.lambda = param.lambda/t;
out = mexProximalTree(x, tree, param);