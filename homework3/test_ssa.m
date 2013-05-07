% --------------------------------
% Test SSA 
%
% Exercise3
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY
% -----------------------------

% load data in X
load('ssa_data.mat');
[d, N] = size(X);

% plot components
%i=10
%plot(X(i,:));
%title(strcat( strcat('x',num2str(i)), '(t)' ))

% define some arguments
Sn_dim = 7;
n_epoches = 25;
epoch_length = N/n_epoches;


% split into epoches defined in epoch_length
epoches = mat2cell( X, [d], ones(1,n_epoches)*epoch_length );

% do ssa
[ Ps, An, y, converged, iter] = ssa_simple(epoches,Sn_dim);

% show projection
Ps

% apply projection Ps to stationary subspace Sn 
Sn =Ps*X;

% plot
plot( Sn(1,:) );
for i=1:Sn_dim
    subplot(Sn_dim,1,i);
    plot(Sn(i,:));
    title(strcat( strcat('s',num2str(i)), '(t)' ))
end