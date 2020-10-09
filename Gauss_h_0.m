function [anh0,numh0] = Gauss_h_0 (sig,n);

MISEh=@(n,h,sig) (-0.5/n - 0.5*(1-1/n)*(1+(sig/h)^2)^(-3/2) + sqrt(2)*(1+(2*sig^2)/h^2)^(-3/2))
epsh =@(n,h,sig) (-2/n - 0.5*(1-1/n)*(1+(sig/h)^2)^(-3/2) + sqrt(2)*(1+(2*sig^2)/h^2)^(-3/2))
nnow=n;
anfun= @(x) MISEh(nnow,x,sig);
numfun=@(x) epsh(nnow,x,sig);
anh0=fzero(anfun,0.5*sig);   % second number is initial guess
numh0=fzero(numfun,0.5*sig); % second number is initial guess

end