function OUT = prox_family_group(x, lambda, C, weights)

if nargin < 4
   weights = ones(size(C,1),1);
end

norms = sqrt(C*x.^2);
norms = max((norms - lambda*weights)./norms,0);
norms = repmat(norms(2:end),1,5)';
norms = norms(:);
norms = [repmat(norms(1),4,1); norms];
OUT = norms.*x;