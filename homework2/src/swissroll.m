% Generate swissroll data
%
% INPUT:
%
%       N (1 x 1): number of points to generate
%
% OUTPUT:
%
%       X (N x 3): coordinates
%       Y (N x 1): true embedding
%
% author(s): Till Rohrmann
function [X, Y] = swissroll(N)
    Y = 2*pi*rand(N,1);
    epsilon = rand(N,1);
    X = [(1+Y).*cos(Y),(1+Y).*sin(Y),epsilon];
end