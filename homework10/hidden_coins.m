function [states,output] = hidden_coins( n )
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY

trans = [.1 0.9;0.5 0.5];
emis = [0.2 0.8; 0.4 0.6];


states = zeros(1,n);
output = zeros(1,n);

% fist default state
states(1)=1;
output(1)=(rand > emis(1,1))+1;

for i=2:n
    x = trans(states(i-1),:);
    next_state=(rand > x(1));
    states(i)=(next_state+1);
    output(i)=(rand > emis(states(i),1))+1;
end
