function sheet02_solution

    % generate a swissroll data set
    if 1 %set to 1 if swissroll.m is done.
      [X, Y] = swissroll(1000);

      figure('Name','swiss roll');
      scatter3(X(:, 1), X(:, 2), X(:, 3), 30, Y, 'filled')

      % apply t-SNE. uncomment to use.
      no_dims = 2;
      init_dims = 3;
      
      
      perplexities = 10:10:90;
      counter = 1;
      m = ceil(sqrt(numel(perplexities)));
      n = ceil(numel(perplexities)/m);
      
      figure();
      for perplexity = perplexities
          mappedX = tsne(X,[],no_dims, init_dims, perplexity);
           % plot mapped data
          subplot(m,n,counter);
          scatter(mappedX(:, 1), mappedX(:, 2), 20, Y);
          title(strcat(['Perplexity:',int2str(perplexity)]));
          counter = counter + 1;
      end
      
      suptitle('mapped swiss roll data');

      %mappedX = X; % use this to show the swiss roll

    end

    % load mnist data
%     load mnist_train.mat
%     rand('twister', 5489)
%     ind = randperm(size(train_X, 1));
%     train_X = train_X(ind(1:1000),:);
%     train_labels = train_labels(ind(1:1000));
% 
%     if exist('mnist_results.mat', 'file')
%       load mnist_results.mat
%     else
%       no_dims = 2;
%       init_dims = 30;
%       perplexity = 30;
% 
%       mappedX = tsne(train_X, [], no_dims, init_dims, perplexity);
% 
%       save mnist_results.mat mappedX
%     end
% 
%     % plot data
%     figure('Name','mnist_embedding');
%     show_mnist(mappedX, train_X, train_labels);

end

