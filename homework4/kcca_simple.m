function [r,a,b] = kccpa_simple( kX, kY, kappa )
% author:  ROHRMANN TILL, FLEISCHMANN KAY
    % dimensions of the data
    n = size(kX,1);
    
    options.disp = 0;

    LH = [zeros(n) kX*kY';kY*kX' zeros(n)];
    RH = [kX*kX' zeros(n);zeros(n) kY*kY'] ;

    % compute the generaized eigenvalue problem
    [Vs,r]=eigs(LH,RH,1,'LA',options);

    % cut off the two different projections
    a = Vs(1:n);
    b = Vs(n+1:end);    
end