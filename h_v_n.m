sig=15;
n=logspace(1,6,61);
clear h

%depsdh=@(n,h,sig) -sqrt(2)/n - 0.5*(1-1/n)*(1+(sig/h)^2)^(-3/2) + sqrt(2)*(1+(2*sig^2)/h^2)^(-3/2))
dmisedh=@(n,h,sig) (-sqrt(2)/n - 0.5*(1-1/n)*(1+(sig/h)^2)^(-3/2) + sqrt(2)*(1+(2*sig^2)/h^2)^(-3/2))
depsdhstab=@(n,h,sig) -1/n - (1-1/n)*(1+sig/h)^-2 + 4*(1+2*sig/h)^-2;

for i=1:length(n)
    nnow=n(i)
    fun= @(x) dmisedh(nnow,x,sig);
    funstab= @(x) depsdhstab(nnow,x,sig);
   
    h(i)=fzero(fun,0.5*sig)   % second number is initial guess
    hstab(i)=fzero(funstab,0.5*sig) ;
end
figure(5)
loglog(n,h/sig,'-d')
hold on
loglog(n,hstab/sig,'-sq')
loglog(n,1.0*n.^(-2./6))
loglog(n,1.06*n.^(-1/5))
axis([10 10^6 .01 2]);

blahh=linspace(0,2,100);
for i=1:length(blahh)
    eps(i)= dmisedh(10^3,blahh(i),sig)  
end
figure(6)
plot(blahh,eps)

[val,idx]=min(eps)