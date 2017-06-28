function OUT = gradOperator_luca(X,m)

X = reshape(X,[m,m]);

P1 = X(1:m-1,:)-X(2:m,:);

P2 =X(:,1:m-1)-X(:,2:m);

OUT = [P1(:); P2(:)];

% scale = [1, 1];
% f = reshape(f,[m,m]);
% [Dx, Dy] = gradOperator(f, scale);
% OUT = [Dx(:); Dy(:)];
% 
% end
% 
% function [Dx, Dy] = gradOperator(f, scale)
% 
%     % Compute the gradient of the images.
%     Dx = scale(1)*[diff(f,1,2), f(:,1) - f(:,end)];
%     Dy = scale(2)*[diff(f,1,1); f(1,:) - f(end,:)];
% 
% end