% Experiment

A = [0.1 0.9;0.5 0.5;];
B = [0.2 0.8;0.4 0.6;];
PI = [1 0];
% observation
y = [1 1 2 2 2 1 1 ];

x = viterbi( A,B, PI, y )
