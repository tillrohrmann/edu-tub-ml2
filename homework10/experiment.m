% Experiment
%
% Authors:  Kay Fleischmann
%           Till Rohrmann
close all;
clear all;
clc;
A = [0.1 0.9;0.5 0.5;];
B = [0.2 0.8;0.4 0.6];
PI = [1 0];

%test viterbi
%observation
%y = [1 1 2 2 1 1 1 2 2 2 1 2 1 1 2 2 1 1 2 2 1 2 1 ];
%x1 = viterbi( A,B, PI, y );
%seq1 = hmmtrain(100,A,B);
%x2 = hmmviterbi(y,A,B);
%sum(x1==x2)/length(y)

% do the experiment
l = [5 10 20];
N = 1000;
for i=1:length(l)
    
    % random
    random_states = rand(N,l(i));
    random_states(random_states>=0.5)=1;
    random_states(random_states<0.5)=2;
    
    stationary_distribution = [0.3571, 0.6429];
    randomNumbers = rand(N,l(i));
    stationary_states = zeros(N,l(i));
    stationary_states(randomNumbers < stationary_distribution(1)) =1;
    stationary_states(randomNumbers >= stationary_distribution(1)) = 2;
    
    viterbi_results = zeros(N,l(i));
    random_results = zeros(N,l(i));
    stationary_results = zeros(N,l(i));
    
    % generate sequences and prove results
    % viterbi
    for n=1:N
        [states,output] = hidden_coins( A,B, l(i) );
        seq = viterbi( A,B, PI, output );

        viterbi_results(n,:) = states == seq;
        random_results(n,:) = states == random_states(n,:);
        stationary_results(n,:) = states == stationary_states(n,:);
        
    end
    
    corr_viterbi = sum(viterbi_results) / N;
    corr_random = sum(random_results) / N;
    corr_stationary = sum(stationary_results) /N;

    total_corr_viterbi = sum( sum(viterbi_results,2)==l(i)) / (N);
    total_random_viterbi = sum( sum(random_results,2)==l(i)) / (N);
    total_stationary_viterbi = sum(sum(stationary_results,2)==l(i))/N;
    
    subplot(length(l),1,i); 
    hold on;
    plot(corr_viterbi)
    
    title(['sequence length l= ' int2str( l(i) ) ', whole sequence viterbi: ' num2str(total_corr_viterbi) ', whole sequence random: ' num2str(total_random_viterbi) ', whole sequence stationary: ' num2str(total_stationary_viterbi)]);
    

    plot(corr_random,'red');
    plot(corr_stationary,'black');
    legend('Viterbi','Random', 'Stationary')

    hold off;
    
end
