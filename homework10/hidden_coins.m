function [states,output] = hidden_coins( trans, emis, n )
%
% Author:   ROHRMANN, TILL
%           FLEISCHMANN, KAY

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
