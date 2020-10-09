function [h0] = FT_h(data)

ndata=length(data);
mu=mean(data);
sig=std(data);

nf=2^18  % # of frequencies
iqrvec=quantile(data,[0.25 0.75]);  %interquartile range of data
iqr=iqrvec(2)-iqrvec(1)
dtF=.1*iqr  % Manual override!

F=zeros(1,nf);

% Make a vector of frequencies
F(1:nf/2+1) = 1*(0:(nf/2))/nf/dtF;  
F(nf/2+2:nf) = -F(nf/2:-1:2);
dF = F(2)-F(1);

% Construct FT based on actual data points (Pankavich) 
for f=0:nf-1
    Fnow=F(f+1);
    Fn2(f+1)=sum(exp(-2*pi*1i*Fnow*(data-mean(data))))/ndata; 
end

GK=@(h,F) exp(-h^2*(2*pi*F).^2/2);  % Get the Gaussian kernel 

nh=1000; h=linspace(0.0,2*sig,nh);
for j=1:nh
    hnow=h(j);
    K=GK(hnow,F);
    eps2(j)=2*dF*sum(K)/ndata + dF*sum( ( (1-1/ndata)*K.^2 - 2*K).*Fn2.^2 );
end
[val,idx]=min(real(eps2));

h0=h(idx)

figure(77)
plot(h,eps2)
hold on
plot(h0,eps2(idx),'+')
h0=min([h0 10*iqr])
xlabel('bandwidth (h)')
ylabel('eps(h)')
title('h_0 at minimum of eps(h)')
end

