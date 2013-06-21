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
Z = zeros(1,lenT); % output

% init
T1 = zeros(numX, lenT ); %store max probability from t-1
T2 = zeros(numX, lenT ); % track best path
t=1; % process first time step
T1(:,t) = pi' .* B( :, y(t));

for t=2:lenT
    for j=1:numX
        [v,k] = max(T1(:,t-1) .* A(:,j) .* B(j,y(t)) );
        T1(j,t) =v;
        T2(j,t) =k;
    end
    
end


[v,k]=max( T1(:,lenT) );
Z(lenT)=k;
x(lenT)=Z(lenT);

% find the best path
for t=lenT:-1:2
    Z(t-1)=T2(Z(t),t);
    x(t-1)=Z(t-1);
end

%T1
%T2

end
