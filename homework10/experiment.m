% Experiment
A = [0.1 0.9;0.5 0.5;];
B = [0.2 0.8;0.4 0.6;];
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
    
    viterbi_results = zeros(N,l(i));
    random_results = zeros(N,l(i));
    
    % generate sequences and prove results
    % viterbi
    for n=1:N
        [states,output] = hidden_coins( l(i) );
        seq = viterbi( A,B, PI, output );

        viterbi_results(n,:) = states == seq;
        random_results(n,:) = states == random_states(n,:);
        
    end
    
    corr_viterbi = sum(viterbi_results) / N;
    corr_random = sum(random_results) / N;

    subplot(length(l),1,i); 
    hold on;
    plot(corr_viterbi)
    
    title(['sequence length l= ' int2str( l(i) )])
    

    plot(corr_random,'red');
    hold off;
    
    legend('Viterbi','Random')
    
end