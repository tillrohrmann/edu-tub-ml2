function [c,U,V,c_perc,U_perc,V_perc] = tkcca(X,Y,tau,kappas)
% [c, U, V, c_bstrp, U_bstrp, V_bstrp] = tkcca(X,Y,tau,kappas)
%
% 	temporal kernel CCA
%	
% INPUT
%	X	a nDimensions-by-Samples data matrix
%	Y	a mDimensions-by-Samples data matrix
%	tau	the maximal time lag by which X is shifted with respect to Y
%	kappas	the regularizers
%	
% OUTPUT
%	c	the canonical correlogram, from -tau to tau
%	U	the time resolved canonical variate for X
%	V	the canonical variate for Y
%	
%	USAGE 
%
%>> [c,U,V] = tkcca(X,Y,10,.1)
%		
%	computes canonical correlogram c, a time resolved variate U and canonical 
%	projection V for timelag -10 to 10 with fixed regularizer of .1
%

if nargin<4, kappas = [10.^-[0:7]]'; end

% center data
X = X - repmat(mean(X,2),1,size(X,2));
Y = Y - repmat(mean(Y,2),1,size(Y,2));
% embed one signal in time shifted copies of itself
[X, timeidx, tauidx] = embed(X,tau);
% compute the linear kernels
kY = Y(:,timeidx)' * Y(:,timeidx);
kX = X' * X;
% find the right regularizer
kappaOpt	= optimize_kappa(kX,kY,kappas);
% compute kcca using the right regularizer
[r, a, b]	= kcca(kX,kY,kappaOpt);
% reconstruct the canonical variates and compute canonical correlogram
[U, V, c]	= reconstruct(X,Y(:,timeidx),a,b,tauidx);


function [eX, timeidx, tauidx] = embed(X,tau)
	% embed the first signal in its temporal context
	[D T] = size(X);
	% in case tau is a scalar, make it a vector from -tau to tau
	if length(tau)==1,tau = -tau:tau;end
	startInd 	= abs(tau(1)) + 1;
	stopInd		= T - abs(tau(end));
	len			= stopInd - startInd + 1;
	% create a column vector that contains the indices of the first segment
	idx = repmat((startInd:stopInd)', 1, length(tau)) + repmat(tau, len, 1); 	
	% create (linear) indices for the different dimensions
	dim_offset = repmat( (0:D-1)*T, length(tau)*len, 1);
	idx = repmat(idx(:), 1, D) + dim_offset;
	% for the linear indices we need column-signals
	X = X';
	% get the data (D channels, segments are concatenated) and reshape it
	eX = reshape(X(idx), len, length(tau)*D)';
	tauidx = repmat(tau',D,1); 
	timeidx = startInd:stopInd;

function [r,a,b] = kcca(kX,kY,kappa)
	% compute the dual coefficients 
	n = size(kX,1);
	options.disp = 0;
	% force kernel symmetry
	kX = (kX+kX')/2; kY = (kY+kY')/2;
	% normalise the spectral norm of the matrices
	kX = kX./max(eig(kX));
	kY = kY./max(eig(kY));
	% Generate LH
	LH = [zeros(n) kX*kY';kY*kX' zeros(n)];
	% generate RH with regularization ridge
	RH = [kX*kX' zeros(n);zeros(n) kY*kY'] + eye(2*n)*kappa;
	% make sure the matrices are symmetric
	RH=(RH+RH')/2; LH=(LH+LH')/2;
	% Compute the generalized eigenvectors
	[Vs,r]=eigs(LH,RH,1,'LA',options);
	a = Vs(1:n);
	b = Vs(n+1:end);
	
function kappa = optimize_kappa(kX,kY,kappas,iterations)
	if nargin<4, iterations=10; end
	shcors = zeros(length(kappas),iterations);
	skX = kX; skY = kY;
	% try each regularizer
	for iR = 1:length(kappas)
		r(iR,1) = kcca(kX,kY,kappas(iR));
		% for all iterations of the reshuffling procedure
		for iS = 1:iterations
			idx = randperm(size(kX,1));
			skX = kX(idx,idx);		
		    % do cca on shuffled data
     		shcors(iR,iS) = kcca(skX,skY,kappas(iR));
		end
	end
  	% pick that regularizer that maximizes the distance between 
  	% true correlations and shuffled data correlations
  	[val,pick]  = max(mean((repmat(r,1,(iterations))-shcors).^2,2));
  	fprintf('Picked kappa=[ %2.8f ]\ncorrelations  %0.2f (true) vs. %0.2f (shuffled)\n',...
  			kappas(pick),r(pick),mean(shcors(pick,:)))
  	kappa = kappas(pick);

function [U,V,c] = reconstruct(eX,Y,a,b,tidx)
	% the number of time lags
	nTau = length(unique(tidx));
	% the number of dimensions
	D = size(eX,1)/nTau;
	% the time-lag sorted indices
	[sorted, sortInds] = sort(tidx);
	% the zero lag canonical component
	pY = ((Y * b)' * Y)';
	% the time shifted canonical components
	pX = repmat(eX * a, 1, size(eX, 2) ) .* eX;
	if D>1
		pX = squeeze(mean(reshape(pX(sortInds,:), [D,nTau, size(eX,2)] )));
	end
	% the correlations between the zero lag component of Y and the 
	% time shifted components of X (i.e. the canonical correlogram)
	c = corr(pY,pX');
	% the time resolved variate
	U = reshape(eX * a, nTau, D)';
	% the other variate
	V = Y * b;
