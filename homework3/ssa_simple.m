function [ Ps, An, y, converged, iter ] = ssa_simple(X, dd, nreps, max_iter, quiet, nomeans)
%SSA_SIMPLE      Stationary Subspace Analysis.
%
%usage 
%  [Ps, An, y, converged, iter] = ssa_simple(X, dd, { nreps: 5,  max_iter: 100, quiet: false, nomeans: false })
%
%input
%  X          Data in one of two possible formats:
%               1. cell array where each X{i} is a (d x n_i)-dataset
%               2. cell array with two elements: 
%                   X{1} is a (d x n) matrix of epoch means and
%                   X{2} is a (d x d x n) array of epoch covariance matrices 

%  dd         Dimensionality of stationary subspace
%  nreps      Optional: number of restarts w/ different init (default: none) 
%  max_iter   Optional: maximum number of iterations (default: 100)
%  quiet      Optional: suppress output (default: false) 
%  nomeans    Optional: Perform SSA w.r.t. to the covariance matrix only. 
%
%output
%  Ps         Projection matrix to stationary sources (dd x d)
%  An         Estimated basis of the non-stationary subspace (d x d)
%  y          Minimum objective function value
%  converged  True, if the optimization has converged
%  iter       Number of iterations
%
%author
%  buenau@cs.tu-berlin.de

% Set default parameters.
if ~exist('max_iter', 'var') || isempty(max_iter), max_iter = 100; end
if ~exist('nreps', 'var'), nreps = []; end
if ~exist('max_iter', 'var'), max_iter = 100; end
if ~exist('quiet', 'var'), quiet = false; end
if ~exist('nomeans', 'var'), nomeans = false; end

% Loop over repetitions and return the solution with the lowest objective function 
% value. 
if ~isempty(nreps)
  r_Ps = cell(1, nreps);
  r_An = cell(1, nreps);
  r_y = zeros(1, nreps);
  r_converged = zeros(1, nreps);
  r_iter = zeros(1, nreps);

  for i=1:nreps
    [r_Ps{i}, r_An{i}, r_y(i), r_converged(i), r_iter(i) ] = ssa_simple(X, dd, [], max_iter, quiet, nomeans);
  end

  [foo, mini] = min(r_y);
  Ps = r_Ps{mini};
  An = r_An{mini};
  y = r_y(mini);
  iter = r_iter(mini);
  converged = r_converged(mini);
  return; 
end

% Sanity check.
if length(X) < 2, error('X must contain at least two datasets or two matrices: X{1} contains the means and X{2} contains the covariance matrices\n'); end

% Distinguish the two parametrization variants: data or means+covariance matrices.
X_contains_data = (ndims(X{2}) == 2);

d = size(X{1}, 1);

% Sanity check.
if dd > d, error('Dimensionality of stationary subspace (dd) must be less than or equal to the dimensionality of the input data!\n'); end

% Distinguish different formats. 

if X_contains_data
  % If the input is data, compute means and covariance matrices per epoch 
  % and the whitening matrix.
  n_X = length(X);
  C = zeros(d, d, n_X);
  mu = zeros(d, n_X);
  all_data = [];
  for i=1:n_X
    all_data = [ all_data X{i} ];
    mu(:,i) = mean(X{i}, 2);
    C(:,:,i) = cov(X{i}'); 
  end
else
  % If the input is means+covmats, compute only the whitening matrix. 
  mu = X{1};
  C = X{2};
  n_X = size(mu, 2);
end

converged = false;
y_new = [];

% Parameters for backtracking linesearch.
ls_alpha = 0.5*(0.01+0.3);
ls_beta = 0.4;

if ~quiet, fprintf('*** iter=[iteration] y=[function value] ||grad||=[norm of gradient] ||step||=[norm of step] ([number of line search iterations]) rel_dec=[relative function value decrease in percent]\n'); end

% Centering and Whitening. 
W = inv(sqrtm(squeeze(mean(C,3))));
mu = mu - repmat(mean(mu,2), [1 n_X]);

% Initialization: random rotation.
B = randrot(d)*W;

% Apply initialization to means and covariance matrices.
mu = B*mu;
C = mult3(C, B);

% Optimization loop.
for iter=1:max_iter
  % Get current function value and gradient.
  [y, grad] = objfun(zeros(dd*(d-dd),1), C, mu, d, dd, nomeans);

  % Sanity check.
  if ~isempty(y_new) && y_new ~= y, error('Something is utterly wrong.\n'); end

  % Print progress (if not suppressed).
  if ~quiet, fprintf('iter=%d y=%.5g ||grad||=%.5g ', iter, y, norm(grad)); end

  % Conjugate gradient: compute search direction.
  if iter == 1
    alpha = -grad;
  else
    gamma = grad'*(grad-grad_old)/(grad_old'*grad_old);
    alpha = -grad + gamma*alpha_old;
  end
  grad_old = grad;
  alpha_old = alpha;

  % Alpha is the current search direction, so we do a line-search along t*alpha
  % Alpha is a vector, that contains the elements of "Z" (see paper)

  % Normalize search direction. 
  alpha = alpha ./ (2*norm(alpha));

  % Fill in nonzero values: this means: put the Z into the bigger M so that we 
  % we can do expm(M)
  M_alpha = reshape(alpha, [dd (d-dd)]);
  M_alpha = [ zeros(dd, dd) M_alpha; -M_alpha' zeros(d-dd, d-dd) ];

  % Backtracking line search loop.
  t = 1;
  for j=1:10
    M_new = t*M_alpha;
    R = expm(M_new);

    y_new = objfun(zeros(dd*(d-dd), 1), mult3(C, R), R*mu, d, dd, nomeans);

    % Stop if function decrease is sufficient.
    if y_new <= (y + ls_alpha*t*grad'*alpha)
      break;
    end

    t = ls_beta*t;
  end

  % Stop if line search failed. 
  if y_new >= y
    if ~quiet, fprintf('no step found\n'); end
    converged = true;
    break;
  end

  % Stop if relative function decrease is below threshold.
  rel_dec = (y-y_new)/y;
  rel_dec_thr = 1e-8;
  if rel_dec < rel_dec_thr
    if ~quiet, fprintf('rel_dec < %f\n', rel_dec_thr); end
    converged = true;
    break;
  end

  % Print progress.
  if ~quiet, fprintf('||step||=%.3g (%d) rel_dec=%.3g%% y=%.3g\n', t, j, 100*rel_dec, y); end

  % Rotate basis (= multiplicative update step).
  C = mult3(C, R);
  mu = R*mu;
  
  B = R*B;
end

% Display warning message if algorithm has not converged.
if ~converged
  if ~quiet, fprintf('Reached maximum number of iterations\n'); end
end

% Compute estimated mixing matrix.
A = inv(B);

% Split estimated de-mixing matrix into two projection matrices.
Ps = B(1:dd,:);
An = A(:,(dd+1):end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = mult3(C, R)
% Compute R*C(:,:,i)*R' for all i.

[d1, d2, d3] = size(C);

% Multiply from the left with R.
C = reshape(C, [d1 d2*d3]);
C = reshape(R*C, [d1 d2 d3]);

% Multiply from the right with R'
C = permute(C, [2 1 3]);
C = reshape(C, [d1 d2*d3]);
C = reshape(R*C, [d1 d2 d3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
function [fx, grad] = objfun(M, C, mu, d, dd, nomeans)

n = size(C, 3);

% Degrees of freedoms.
dof = n*(dd*(dd+1)/2 + dd);

% Fill in non-zero values. 
M = reshape(M, [dd (d-dd)]);
M = [ zeros(dd, dd) M; -M' zeros(d-dd, d-dd) ];

% Compute rotation.
R = expm(M);

% Projection to stationary signals.
P = R(1:dd,:);

log_det = zeros(1, n);
inv_pC = zeros(dd, dd, n);
pC = zeros(dd, dd, n);
log_det = zeros(1, n);
pmu = zeros(dd, n);

opts.UT = true;

% Prepare some values.
for i=1:n
  pC(:,:,i) = P*C(:,:,i)*P';
  L = chol(pC(:,:,i));
  inv_L = linsolve(L, eye(dd), opts);
  inv_pC(:,:,i) = inv_L*inv_L';
  log_det(i) = 2*sum(log(diag(L)));
  pmu(:,i) = P*mu(:,i);
end

fx = 0; 
grad = zeros(dd,d);

% Loop over epochs.
for i=1:n
  fx = fx - log_det(i);
  if ~nomeans, fx = fx + pmu(:,i)'*pmu(:,i); end

  grad = grad -2*inv_pC(:,:,i)*P*C(:,:,i);
  if ~nomeans, grad = grad + 2*P*mu(:,i)*mu(:,i)'; end
end

grad = [ grad; zeros(d-dd, d) ];
grad = grad*R' - R*grad';

grad = grad(1:dd,dd+1:end);
grad = grad(:);

%fx = (fx-dof)/sqrt(2*dof);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ R, M ] = randrot(d)
%RANDROT        Generate random orthogonal matrix. 
%
%usage
%  [R,M] = randrot(d)
%
%author
%  buenau@cs.tu-berlin.de

M = 10*(rand(d,d)-0.5);
M = 0.5*(M-M');
R = expm(M);
