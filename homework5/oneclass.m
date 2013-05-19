function alpha = oneclass(K, C)
% inputs:  K      N x N kernel matrix
%          C      regularization constant
% outputs: alpha  N-dimensional dual solution vector

dim=length(K);

%variables for quadratic program
%minimize   c' * x + 1/2 x' * H * x
%subject to A'*x = b
%           l <= x <= u

%u=ones(dim,1);
%c=diag(K)';
%H=-K;
%b=1;
%A=u;
%l=zeros(dim,1);
%u=ones(dim,1)*C;
%[x,y] = pr_loqo2(c, H, A, b, l, u)

% 1/2*x'*H*x + f'*x 
% A*x ≤ b
% [ ] *x ≤ [C,..,C]'
% Aeq*x = beq
H=-K;
f=diag(K);
l=zeros(dim,1);
u=ones(dim,1)*C;
Aeq = ones(dim,1)';
beq = [1];
alpha = quadprog(H,f,[],[],Aeq,beq,l,u)