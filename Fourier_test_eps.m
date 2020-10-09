clear all

for nens=1:5  % Do an ensemble
ndata=1000;
nbins=2^nextpow2(ndata/10);   % Also defines spatial discret. of histogram-based FT
sig=1; mu=5;
data = -log(rand(1,ndata));  % Exponential
%data = mu+sig*randn(1,ndata);  % Gaussian

%Construct a histogram ... 
% a sum of delta functions mapped to nearest coords.
figure(1)
h=histogram(data,nbins);
dxhist=h.BinWidth
minx=h.BinLimits(1); maxx=h.BinLimits(2); 
X=h.Values/length(data)/dxhist;

dt=16*dxhist;  % This is the real-space discret. for Pankavich's FT.

n=2^nextpow2(length(X));

nf=2^16;   % arbitrary number of frequencies for Pankavich's method
dt = 8*dxhist;
dt=.1*sig  % Manual override!

F=zeros(1,n);

% Make a vector of frequencies
F(1:nf/2+1) = 1*(0:(nf/2))/nf/dt;  
F(nf/2+2:nf) = -F(nf/2:-1:2);
dF = F(2)-F(1);

% Construct FT based on histogram
for f=0:n-1
    expj=exp(-f*2*pi*1i*(1/n)*(0:n-1));
    Ymanual(f+1)=sum(X.*expj);
end

% Construct FT based on actual data points (Pankavich) 
for f=0:nf-1
    Fnow=F(f+1);
    Ymanual2(f+1)=sum(exp(-2*pi*1i*Fnow*(data-mean(data))))/ndata; 
end
Fn=Ymanual;
Fn2=Ymanual2;
% For a Gaussian kernel K, get K(freq), either by fft or direct
GK=@(h,F) exp(-h^2*(2*pi*F).^2/2);
KK=Fn2;   %  The experimantal density IS the kernel (need to work on scale)
nh=2000; h=linspace(0.0,3*sig,nh);
for j=1:nh
    hnow=h(j);
    K=GK(hnow,F);
%    eps(j)=2*dF*sum(K)/ndata + dF*sum( ( (1-1/ndata)*K.^2 - 2*K).*Fn.^2 );    
    eps2(j)=2*dF*sum(K)/ndata + dF*sum( ( (1-1/ndata)*K.^2 - 2*K).*Fn2.^2 );
    epsK(j)=2*dF*sum(KK)/ndata + dF*sum( ( (1-1/ndata)*KK.^2 - 2*KK).*Fn2.^2 );

end


% Find and plot h=argmin(eps(h))
[val,idx]=min(real(eps2));

hens(nens)=h(idx);

figure(4)
plot(h,real(eps2),'k-')
hold on
%plot(h,real(epsK),'o');
plot(h(idx),val,'b+')


end   % End of ensemble loop

epsGauss=(2/ndata/sqrt(2*pi))./h + ...
1/2./sqrt(pi*(h.^2+sig^2)) - ...
1./sqrt(pi*(0.5*h.^2+sig^2)) - ...
1/2/ndata./sqrt(pi*(h.^2+sig^2));
[valGauss,idxGauss]=min(real(epsGauss));

taylor=1.06*sig*ndata^-0.2
exact=h(idxGauss)
numerical=mean(hens)
numsd=std(hens)

figure(4)
plot(h,epsGauss,'r');
plot(h(idxGauss),valGauss,'r+')
%axis([0 3*sig -0.4 0])
plot([numerical-numsd numerical+numsd],[0 0],'-+')
plot(taylor,0,'o')



%Ymanual
%Y = fft(X)
%Ymanual2
Xmanual2=ifft(Ymanual2)/dt;
Xmanual=ifft(Ymanual);

% Not sure why I have to shift the Pankavich one but not histogram below
figure(2)
plot(dxhist*(0.5+(0:n-1))+minx,real(Xmanual),'-o')
axis([2*minx-10 2*maxx min(real(Xmanual2)) 1.2*max(real(Xmanual)) ])
hold on
plot(data,zeros(size(data)),'o')
plot(dt*(-nf/2:nf/2-1),real(fftshift(Xmanual2)),'k-+')
hold off

