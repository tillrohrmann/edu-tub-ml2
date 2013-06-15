% Experiment Hidden Markov Model
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY


% Generate 1000 times a random state-sequence and observations
N=1000;
len=20;
X=zeros(N,1);
tails=0;
heads=0;
for i=1:N
    [states, output] = hidden_coins(len);
    tails=tails+sum(output==2)/len;
    heads=heads+sum(output==1)/len;
    X(i)=heads/tails;
end
plot(X);