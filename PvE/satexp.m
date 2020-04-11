function [ Ypred ] = satexp( beta, X )
%saturating exponential with y-offset

Xk   = beta(1);
Max  = beta(2);

Ypred = Max.*(1-exp(-X./Xk));

end

