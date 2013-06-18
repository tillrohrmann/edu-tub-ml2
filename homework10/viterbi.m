function x = viterbi(A, B, pi, y)
% input     A   tansitions between states
%           B   emission probabilties
%               from each state
%           pi  initial state probabilities
%           y   observed sequences
% output    most likely sequences
% 
% Author: KAY FLEISCHMANN,TILL ROHRMANN
%
numX=length(A); % num states
numY=size(B,2); % num possible observations
lenT=length(y); % time
x = zeros(1,lenT); % output

% init
T1 = zeros(numX, lenT ); %store max probability from t-1
T2 = zeros(numX, lenT ); % track best path
t=1; % process first time step
T1(:,t) = pi' .* B( :, y(t));

for t=2:lenT
    [val,index] = max ( (repmat( T1(:,t-1)', numX, 1 )  .* A) , [], 2) ;
    T1(:,t) = val  .* B( :, y(t) );
    T2(:,t) = index;
end
%T1
%T2
[v,z]=max( T1(:,lenT) );
x(lenT)=T2(z,lenT);

% find the best path
for t=lenT:-1:2
    [v,z]=max( T2(z,t) );
    
    x(t-1)=T2(z,t);
end

%T1
%T2

end
