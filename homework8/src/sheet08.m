function sheet08 

X = textread('splice-train-data.txt', '%s');
Y = load('splice-train-label.txt');
XE = textread('splice-test-data.txt', '%s');
YE = load('splice-test-label.txt');

tic
fprintf('Writing training features for k = 1\n');
write_wdk_features('wdk1-train.txt', 1, X, Y);
fprintf('Writing test features for k = 1\n');
write_wdk_features('wdk1-test.txt', 1, XE, YE);
toc

tic
fprintf('Writing training features for k = 2\n');
write_wdk_features('wdk2-train.txt', 2, X, Y);
fprintf('Writing test features for k = 2\n');
write_wdk_features('wdk2-test.txt', 2, XE, YE);
toc

tic
fprintf('Writing training features for k = 3\n');
write_wdk_features('wdk3-train.txt', 3, X, Y);
fprintf('Writing test features for k = 3\n');
write_wdk_features('wdk3-test.txt', 3, XE, YE);
toc

end

function beta = beta(K, k)
beta = 2 * (K - k + 1) / K / (K + 1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Your solution below
%
% Authors: Kay Fleischmann
%          Till Rohrmann

% 1. write out weighted degree kernel features
% out into file FN up to degree K. (You should accept values of K = 1, 2,
% 3
function write_wdk_features(FN, K, X, Y)
    A = {'A','C','G','T'}';
%     A={'A','C'}';
    temp={''};
    X=cell2mat(X);
%     X=X(1:2,1:10);
    n=size(X,2);
    l=size(X,1);
    
    result = [];
    
    for i=1:K
        temp = allcomb(temp,A);
        indices = reshape((repmat(1:i,n-i+1,1)+repmat((0:(n-i))',1,i))',1,(n-i+1)*i);
        substrings = reshape(X(:,indices)',[i,(n-i+1),l]);
        pattern = cell2mat(temp)';
        
        substrings = reshape(repmat(substrings,1,length(temp)),[i,(n-i+1),length(temp),l]);
        pattern = reshape(repmat(pattern,(n-i+1),l),[i,(n-i+1),length(temp),l]);
        
        matches = reshape(sum(substrings==pattern,1)==i,(n-i+1)*length(temp),l);
        
                
        beta = sqrt(2*(K-i+1)/(K*(K+1)));
        
        result = [result,beta*matches'];
    end
    
    [r,c,v] = find(result');
    
    fid = fopen(FN,'w');
    
    if fid == -1
        error('Could not open file');
    end
    
    for i=1:l
        fprintf(fid,'%d',Y(i));
        fprintf(fid,' %d:%f',[r(c==i)';v(c==i)']);
        fprintf(fid,'\n');
    end
end

function C=allcomb(A,B)
    indicesA = 1:length(A);
    indicesB = 1:length(B);
    
    [combinedA,combinedB] = ndgrid(indicesA,indicesB);
    
    C = cellstr(cell2mat([A(combinedA(:)),B(combinedB(:))]));
end