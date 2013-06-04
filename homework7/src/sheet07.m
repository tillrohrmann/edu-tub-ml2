function sheet07

%
% SINC DATA
%
figure(1)

% generate some data
[X, Y] = sincdata(100, 0.1);

% generate the kernel matrix
K = rbfkern(0.1, X);
 
% Run rde
[D, Yh] = rde(K, Y);

% plot the data
plot_fit(X, Y, Yh);
title(sprintf('sinc data set, effective dimensionality = %d', D));

%
% SINE DATA
%
figure(2)

% generate some data
[X, Y] = sinedata(100, 4);

% generate the kernel matrix
K = rbfkern(0.1, X);
 
% Run rde
[D, Yh] = rde(K, Y);

% plot the data
plot_fit(X, Y, Yh);
title(sprintf('sine data set, effective dimensionality = %d', D));

end

function K = rbfkern(w, X)
N = size(X, 1);
XX = sum(X.*X, 2);
D = repmat(XX, 1, N) + repmat(XX', N, 1) - 2 * X * X';
K = exp(-D/(2*w));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Insert your solutions below
%
%
% Authors: 
%   Kay Fleischmann, 352247
%   Till Rohrmann, 343756

function [X, Y] = sincdata(N, r)
    a = -4;
    b = 4;
    X = rand(N,1)*(b-a)+a;
    Y = sin(pi*X)./(pi*X) + r*randn(N,1);
end

function [X, Y] = sinedata(N, K)
    a = -pi;
    b = pi;
    X = rand(N,1)*(b-a)+a;
    Y = sin(K*X)+0.3*randn(N,1);
end

function [D, Yh] = rde(K, Y)
    n = size(K,1);
    drange = (1:floor(n/2.))';
    
    [eigenvectors, eigenvalues] = eig(K);
    eigenvalues = diag(eigenvalues);
    
    % sort eigenvectors in descending order
    [~,idx] = sort(-eigenvalues);
    eigenvectors = eigenvectors(:,idx);
    
    S = eigenvectors'*Y;
    SS = S.*S;
    pSS = cumsum(SS);
    pSS = pSS(1:floor(n/2.));
    rpSS = cumsum(flipud(SS));
    rpSS = flipud(rpSS(floor(n/2.):n-1));
    ML=(drange/n).*log(pSS./drange) + (n-drange)/n.*log(1./(n-drange).*rpSS);
    
    [~,D] = min(ML);
    
    Yh = eigenvectors(:,1:D)*(eigenvectors(:,1:D)'*Y);
end

function plot_fit(X, Y, Yh)
    hold on;
    scatter(X,Y,40,'fill');
    scatter(X,Yh,40,'r','fill');
    legend('Input data','Denoised data');
    hold off;
end