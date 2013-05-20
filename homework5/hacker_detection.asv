load('stud-data.mat')
 
 % compute kernel matrices
 disp('computing kernel matrices...')
 KR = full(Xtr'*Xtr);
 KS = full(Xts'*Xts);
 KSR = full(Xts'*Xtr);
 
 % compute the alphas
 disp('learning one-class-SVM...')
 C = 0.002; % adjust C 
 alpha = oneclass(KR, C);
 
 % compute anomaly scores
 as = compute_scores(KS, KSR, KR, alpha); 
 
 Ap = (as > 1);
 
 predicted_attacks = find(Ap)';
 
 hold on;
 plot(as);
 plot(1:length(as),1,'r-');
 
 predicted_attacks