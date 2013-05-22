function scores = compute_scores(KS, KSR, KR, alpha)
% inputs: KR      kernel matrix on training data
%         KS      kernel matrix on test data
%         KSR     kernel matrix on test data/training data
%         alpha   learnt dual vector
% output: scores  vector of outlier scores
N=length(KS);
scores = zeros(N,1);
for i=1:N
    scores(i,1)=KS(i,i)-2*KSR(i,:)*alpha+alpha'*KR*alpha;
end