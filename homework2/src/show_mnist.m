% Show results on the mnist data set. Show all points as dots,
% show 50 random samples as small image.
%
% INPUT:
%
%       mappedX         (N x d): the data mapped to R^d
%       train_X         (N x 784): pixel vectors of 28 x 28 images
%       train_labels    (N x 1): class labels
%
% author(s): Till Rohrmann
function show_mnist(mappedX, train_X, train_labels)
    n = 100;
    indices = randperm(size(train_X,1));
    indices = indices(1:n);
    imSize = 0.04;
    
    hold on;
    h =scatter(mappedX(:,1),mappedX(:,2),23,train_labels);
    %set(gca,'XTickLabel','', 'YTickLabel','');
    
%     xvalues = xlim;
%     yvalues = ylim;
%     position = get(gca,'position');
%     
%     c = colormap();
%     colormap([c; flipdim(gray(64),1)]);
%     
%     cmin = min(train_labels);
%     cmax = max(train_labels);
%     m = 64;
%     C = min(m,round((m-1)*(train_labels-cmin)/(cmax-cmin))+1);
%     
%     set(h,'CData',C);
%     caxis([1,128]);
%     for index = indices
%        image = reshape(train_X(index,:),28,28);
%        axes('position',[position(3)*(mappedX(index,1)-xvalues(1))/(xvalues(2)-xvalues(1)) + position(1)-0.5*position(3)*imSize...
%            position(4)*(mappedX(index,2)-yvalues(1))/(yvalues(2)-yvalues(1))+position(2)-0.5*position(4)*imSize...
%            ,position(3)*imSize,position(4)*imSize])
%        image = flipdim(image',1);
%        h = pcolor(image);
%        C = min(m,round((m-1)*image)+1);
%        C = C +64;
%        set(h,'CData',C);
%        shading flat;
%        axis off;
%        caxis([1,128]);
%     end
    
    hold off;
end