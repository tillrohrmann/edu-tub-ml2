% Experiment Hidden Markov Model
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY


% Generate 1000 times a random state-sequence and observations
N=1000;
len=20;
tails=zeros(len,N);
heads=zeros(len,N);
for i=1:N
    [~, output] = hidden_coins(len);
    tails(:,i)=output==2;
    heads(:,i)=output==1;
end

means = mean(heads,2);
stds = std(heads,0,2)/sqrt(N);

errorbar(1:len,means,stds,'o');
title('Relative freq. to observe heads plus its standard deviation');
xlabel('n');
ylabel('Relative freq.');