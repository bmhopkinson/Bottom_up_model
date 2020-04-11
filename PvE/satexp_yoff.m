function [ Ypred ] = satexp_yoff( beta, X )
%saturating exponential with y-offset

Xk   = beta(1);
Max  = beta(2);
Yoff = beta(3);

Ypred = Max.*(1-exp(-X./Xk)) + Yoff;

end

