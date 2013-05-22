%load('stud-data.mat')
 
% compute kernel matrices
%disp('computing kernel matrices...')
%KR = full(Xtr'*Xtr);
%KS = full(Xts'*Xts);
%KSR = full(Xts'*Xtr);
 
% compute the alphas
%disp('learning one-class-SVM...')

%C=0.002; %adjust C

% find appropriate C
%for i=1:10
%    C=C/2; 
%    alpha = oneclass(KR, C);
%
%    % compute anomaly scores
%    as = compute_scores(KS, KSR, KR, alpha); 
%
%   Ap = (as > 1);
%
%    predicted_attacks = find(Ap)
%    results(i)=length(predicted_attacks);
%   plot(results);
%end

% find hackers in trainingsdata
%alpha = oneclass(KR, C);

% compute anomaly scores
%as = compute_scores(KS, KSR, KR, alpha); 

%Ap = (as > 1);
%predicted_attacks = find(Ap)

% compute anomaly scores
as_t = compute_scores(KR, KSR', KS, alpha); 

Ap_t = (as_t > 1);
predicted_attacks = find(Ap_t)
   

% find hackers in test-data

% compute anomaly scores
%as = compute_scores(KS, KSR, KR, alpha); 

%Ap = (as > 1);
%predicted_attacks = find(Ap)
   

