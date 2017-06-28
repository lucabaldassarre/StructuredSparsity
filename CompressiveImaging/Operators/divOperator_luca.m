function OUT = divOperator_luca(P,m)

m2 = m-1;
n2 = m;
m1 = m;
n1  = m-1;

P1 = P(1:m2*n2);
P1 = reshape(P1,[m2,n2]);
P2 = P(m2*n2+1:end);
P2 = reshape(P2,[m1,n1]);

m = m2+1;
n = n2;

X           = zeros(m,n);
X(1:m-1,:)  = P1;
X(:,1:n-1)  = X(:,1:n-1)+P2;
X(2:m,:)    = X(2:m,:)-P1;
X(:,2:n)    = X(:,2:n)-P2;

OUT = X(:);

% P{1} = 
% scale = [1, 1];
% y = x(m^2+1:end);
% x = x(1:m^2);
% 
% OUT = divOperator(x, y, scale);
% 
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% FUNCTION: DTxyz = divOperator(x, y, scale)
% %%% PURPOSE:  Define the gradient operator of a 2D image f.
% %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DTxy = divOperator(x, y, scale)
% 
%     % Compute the divergence operator of x and y.
%     DTxy = [x(:,end) - x(:,1), -diff(x,1,2)];
%     DTxy = scale(1)*DTxy + scale(2)*[y(end,:) - y(1,:); -diff(y,1,1)];
%     
% end