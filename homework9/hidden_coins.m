function [states,prob] = hidden_coins( n )
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY


trans = [.1 0.9;0.5 0.5];
emis = [0.2 0.8; 0.4 0.6];
initial = [1 0];
states = [1 1; 1 2; 2 1; 2 2;];
states = [1 1 1; 1 1 2; 1 2 1;1 2 2;2 1 1; 2 2 1; 2 2 2; 2 1 2;];
N=length(states);
len=size(states,2);
seq = ones(N,len)*2;

% compute p(x) sequence-transition probabilties
px_seq = horzcat( initial( states(:,1))', trans( sub2ind( size(trans), states(:,1:end-1),  states(:,2:end)) ));
px = prod( px_seq,2);

%compute p(y|x)
pyx = prod( px_seq .* emis( sub2ind( size(emis), seq, states )), 2);

py = sum(pyx.*px);

pxy = (pyx.*px)./py;

prob=pxy;

end