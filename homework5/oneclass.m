function alpha = oneclass(K, C)
% inputs:  K      N x N kernel matrix
%          C      regularization constant
% outputs: alpha  N-dimensional dual solution vector

dim=length(K);

%variables for quadratic program
%minimize   c' * x + 1/2 x' * H * x
%subject to A'*x = b
%           l <= x <= u
% Dimensions: c : N-column vector
%             H : NxN matrix
%             A : N-row vector
%             b : real number
%             l : N-column vector
%             u : N-column vector
% 
%             x : N-column vector
%             y : Objective value
%ERROR: chol works not propperly on very slow (close to zero) data
%c=diag(K)';
%H=-K;
%b=1;
%A=ones(dim,1);
%l=zeros(dim,1);
%u=ones(dim,1)*C;
%[x,y] = pr_loqo2(c, H, A, b, l, u);
%alpha = x;

H=-K;
f=diag(K);
l=zeros(dim,1);
u=ones(dim,1)*C;
Aeq = ones(dim,1)';
beq = ones(1);
alpha = quadprog(H,f,[],[],Aeq,beq,l,u);
