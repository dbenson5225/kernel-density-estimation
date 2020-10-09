clear all
disp('Beginning smoother...');

% Import times and densities from PT simulation
% fid   = load('well.extraction.txt');
% Define some variables
alpha   = 0.5;     % Adaptive factor exponent
Galpha  = 0.5;     %    " for Gauss kernel
ninterp = 12000;    % Number of linearly spaced points to get density estimate
imethod = 'linear' % Matlab's interpolation scheme (linear or spline)
nens    = 1;      % Number in ensemble (choose 1 for BTC data)
itermax = 200;     % Max number of iterations
tol     = 1e-9;    % Inter-iteration density RMSE tolerance 
ndata   = 1000;    % Number of "fake" data to generate
mu=100; sig=10;       % mean and scale of data to be generated

% Placeholders
dkl2=0; Gl2=0; G1passl2=0; Gh0l2=0;
dkl2par=0; Gl2par=0; G1passl2par=0; Gh0l2par=0;

for l=1:nens    % Ensemble loop

%  Generate (or read) some data after choosing type 
%rv='IG';
%rv='Gauss';
%rv='Cauchy';
%rv='Stable';
%rv='Exponential';
%rv='Levy';
rv='readfile'   % Choose this to read a file of data, like BTC data.
if(strcmpi(rv,'readfile'))
      ranvec=load('BTC_50000');
      ranvec=ranvec(ranvec(:,1)>1e-20,1);
      ndata=length(ranvec);
end
if(strcmpi(rv,'Gauss'))  % generate Gaussians
      density =@(mu,sig,tvec) exp(-(tvec-mu).^2./(2*sig^2))/sqrt(2*pi*sig^2);
      ranvec=sig*randn(ndata,1)+mu;  % Generate Gaussians
end
if(strcmpi(rv,'Cauchy'))  % generate Cauchys
      density =@(mu,sig,tvec) (sig/pi)./ ((tvec-mu).^2+sig^2);
      ranvec=mu+sig*tan(pi*(rand(ndata,1)-0.5));  % Generate Cauchys
end
if(strcmpi(rv,'Exponential'))  % generate "shifted" Exponentials
      ranvec=-log(rand(ndata,1))*sig+(mu-sig);  
      %density = @expdens(mu,sig,tvec);
      density = @expdens;
end
if(strcmpi(rv,'Levy'))   % Generate max-skewed gamma-stable RVs
      gamma=0.5; beta=1;
      blah=zeros(ndata,1);
      pd5 = makedist('Stable','alpha',gamma,'beta',beta,'gam',sig,'delta',mu);
      density1 =@(shift,size,tvec) sqrt(sig/2/pi)*exp(-0.5*sig./(tvec-shift))./((tvec-shift).^1.5)
      density = @(shift,size,tvec) pdf(pd5,tvec);
      ranvec=random('Stable',gamma,beta,sig,mu,size(blah));  % Generate 
end
if(strcmpi(rv,'Stable'))   % Generate max-skewed gamma-stable RVs
      gamma=1.5; beta=1;
      blah=zeros(ndata,1);
      pd5 = makedist('Stable','alpha',gamma,'beta',beta,'gam',sig,'delta',mu);
      density = @(shift,size,tvec) pdf(pd5,tvec);
      ranvec=random('Stable',gamma,beta,sig,mu,size(blah));  % Generate 
end
if(strcmpi(rv,'IG'))   % Generate Inverse Gaussian (which may be heavy-tailed)
      density =@(mu,sig,tvec) sqrt(sig./(2*pi*tvec.^3)).*exp(-sig*(tvec-mu).^2./(2*mu^2*tvec));
      yvec=randn(ndata,1);
      yvec=yvec.^2
      xvec=mu + mu^2*yvec/(2*sig) - (mu/2/sig)*sqrt(4*mu*sig*yvec + mu^2*yvec.^2)
      uni=rand(ndata,1);
      apple=find(uni>mu./(mu+xvec));
      xvec(apple) = mu^2./xvec(apple);
      ranvec=xvec;
end

fid(1:ndata,3)=1;      % This is a structure useful to importing particle arrival times
fid(1:ndata,2)=ranvec; % So I'll keep it
fid=sortrows(fid,2);   % Need to sort data so that the kernel can be interpolated later

tvec  = fid(:,2);       % X-axis values
npart = fid(:,3);       % Particle concentration (ones and zeros for observed times/padding)
ndata = size(tvec,1);   % Actual number of data 

% Generate an interpolation grid based on min and max of data ...
texpand=1; tspan=max(tvec)-min(tvec);  % Need to expand range to make kernel extrapolate
tgrid=linspace(min(tvec)-(texpand*tspan),max(tvec)+1*(texpand*tspan),ninterp);  % regularly spaced points to calc. f(t).

% Or manually generate interpolation grid (might be wanted for BTCs)
if(strcmpi(rv,'readfile'))
tgrid=logspace(0,6,ninterp);
%tgrid=linspace(0,1e6,ninterp);
end
dt=diff(tgrid);

%****************
% Generate arrays
%****************

lambda   = zeros(ndata,1);     
h_2t     = zeros(ndata,1);
p_0t    = zeros(ninterp,1);
p_1t    = zeros(ninterp,1);      
p0      = zeros(ninterp,1);
p1      = zeros(ninterp,1);
p1new   = zeros(ninterp,1);
kernel  = zeros(ninterp,1);

%****************
% Solve variables
%****************

%  Get an initial estimate of mean and variance via trapezoidal integration for a histogram ...

%mass     = 0.5*sum(npart(1:end-1)+npart(2:end).*diff(tvec))
%meantimes= 0.5*(tvec(1:end-1)+tvec(2:end));
%center   = 0.5*sum((npart(1:end-1)+npart(2:end)).*diff(tvec).*meantimes)/mass
%variance = 0.5*sum((npart(1:end-1)+npart(2:end)).*diff(tvec).*(meantimes-center).^2)/mass
%width    = sqrt(variance)/1
%scalemin=hmin/width;

% ... or raw data not binned (i.e., raw particle arrival times)
apple=find(npart>0);
mass=sum(npart);
center=sum(npart(apple).*tvec(apple))/mass
datavariance=sum((tvec(apple)-center).^2)/mass
iqr=1*(tvec(apple(floor(0.75*ndata))) - tvec(apple(floor(0.25*ndata))))  %interquartile range
width=iqr/1.35;
width=iqr/1.1;  % Exponential
width=iqr/2.0;  % Cauchy
width=iqr/1.5  % Sort of average?
%width=sqrt(variance)
%scalemin=hmin/width

% Use these for a first approximation:  The width of the data will be the
% "concrete" width against which the kernel is scaled.  Choose one method:
%h_0=0.1*width
h_0 = FT_h(tvec)               % 1) Use Fourier transform routine
%h_0=min([h_0 width])
%h_0 = 1*sig*1.06*ndata^-0.2    % 2) Taylor series approx.
%[h0an,h0num] = Gauss_h_0 (sig,ndata);
%h_0 = h0num     % 3) Exact for Gaussian kernel and a realiz. of Gauss data 
%h_0 = h0an      % 4) Exact for Gauss kernel AND exactly Gauss data

Gh_0=h_0;       % Use a separate h_0 for the Gaussian kernel

%%%%%%%%%%%%%%%%%%%%%%%   Now make some density interpolations %%%%%%%%%%%%%%%%%%%%%

%  First try, using Gaussian kernel ...
parfor j = 1:ninterp
    tnow=tgrid(j);
    p_0t(j) =(1/mass)*sum( npart./sqrt(2*pi*(h_0)^2) .* exp( -(tvec-tnow).^2 / (2*(h_0)^2)) );
end

% The adaptive scheme now uses to  density estimate to get h at each data point 
% Need the value of the interpolated density at each data point - might 
% need interpolate from p_0t if grid not equal to data points

p_atdata=interp1(tgrid,p_0t,tvec,imethod,0);
p_atdata(p_atdata<1e-20)=1e-20;    % If spline is used, can be negative density values

%  Find the adaptive kernel bandwidth at each data point based on the Gaussian convolution:
G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-Galpha;
h_1t = max(width/1000, min(1000*width, h_0.*lambda) );  %may need to put real upper bound

% Standard Gaussian kernel to start, using same notation as when I switch to kernel
% based on p(t):

kernel=(1/sqrt(2*pi*1^2))*exp((tgrid-center).^2/(-2*1^2));  

parfor j=1:ninterp    % size of my interpolation grid
    tnow=tgrid(j);
    timesneeded=center+(tnow-tvec)./h_1t.^1.0;     %scaled times from current grid point to each data point
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0);
    p_1t(j)=(1/mass)*sum( npart.*(1./h_1t.^1.0).*kernelinterp );
    %p_1t(j) = sum( npart./sqrt(2*pi*h_1t.^2) .* exp(-(tvec-tnow).^2 ./ (2*h_1t.^2)) );  
end

p_Gauss=p_1t;  % without a better estimate of h, this is the best I can do.

%DB: Assume that p_1t is a good first estimate of the BTC, using Gaussian kernels.
% Now use the BTC itself as the kernel.  

% start with bandwidth ~ Gaussian optimal

mass=sum(npart);
center=sum(npart.*tvec)/mass;
datavariance=sum((tvec(npart>0)-center).^2)/mass
width=sqrt(datavariance);
%  Or try interquartile or similar?
%apple=find(npart>0.5); numdata=length(apple);
iqrvec=quantile(tvec,[0.25 0.75]);
iqr=iqrvec(2)-iqrvec(1);
iqr=1.*(tvec(floor(0.75*ndata)) - tvec(floor(0.25*ndata)));  %interquartile range of data
%width=iqr/1.35  % Gaussian
width=iqr/1.1;   % Exponential
width=iqr/2.0;   % Cauchy
width=iqr/1.5   % Sort of average?
%scalemin=hmin/width;

kernel=p_1t;  % use the first estimate of BTC as the kernel
intdens=trapz(tgrid,kernel);
%intdens=0.5*sum(dt.*(kernel(1:end-1)+kernel(2:end))')
%intdens=sum(kernel)*dt   % check density sums to one

Gkernel=(1/sqrt(2*pi))*exp((tgrid-center).^2/(-2));
intdens=trapz(tgrid,Gkernel)
%intdens=0.5*sum(dt.*(Gkernel(1:end-1)+Gkernel(2:end)))
%intdens=sum(Gkernel)*dt

figure(5)
plot(center+(tgrid-center).*(1/width),(width).*kernel);
hold on;
plot(tgrid,kernel);
plot(tgrid,Gkernel);
legend('"Standard" kernel','Unscaled kernel','Standard Gaussian')
%pause

parfor j = 1:ninterp
    tnow=tgrid(j);
    timesneeded=center+(tnow-tvec).*width./h_0; 
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0) ;
    p1(j)=(1/mass)*sum( npart.*(width./h_0^1.0).*kernelinterp );
end
parfor j = 1:ninterp
    tnow=tgrid(j);
%    timesneeded=center+(tnow-tvec)./h_0.^1.0; 
    timesneeded=center+(tnow-tvec)./h_0.^1.0; 
    Gkernelinterp=interp1(tgrid,Gkernel,timesneeded,imethod,0) ;
    Gp1(j)=(1/mass)*sum( npart.*(1./h_0.^1.0).*Gkernelinterp );
end

%   Iterate kernel shape and size ........................................

rmselast=1e6; iter=0; fail=0;
intdens=trapz(tgrid,p1)
%intdens=0.5*sum(dt.*(kernel(1:end-1)+kernel(2:end))')
kernel=max(0,p1)./trapz(tgrid,p1);   % Normalize kernel to unit area

while rmselast>tol&iter<itermax
%  recalculate the kernel scaling factor

m_0=trapz(tgrid,p1)
%m_0 = 0.5*sum((kernel(1:end-1)+kernel(2:end))'.*dt)
cum=cumsum(0.5*(dt.*(kernel(1:end-1)+kernel(2:end))'));
%blah=tgrid.*kernel';
%m_1 = 0.5*sum((blah(1:end-1)+blah(2:end)).*dt)/m_0
m_1=trapz(tgrid,(tgrid.*kernel'))
%blah=(tgrid-m_1).^2.*kernel';
%m_2 = 0.5*sum((blah(1:end-1)+blah(2:end)).*dt)/m_0
m_2=trapz(tgrid, (tgrid-m_1).^2.*kernel'  ) 

q3=min(find(cum>0.75));
q1=max(find(cum<0.25));
iqr=(tgrid(q3)-tgrid(q1));
%width=sqrt(m_2) %+ std(tvec)  % width of the kernel
width=iqr/1.35;   % Gaussian
width=iqr/1.1;    % Exponential
width=iqr/2.0;    % Cauchy
width=iqr/1.5    % Sort of average?
iter=iter+1

p_atdata=interp1(tgrid,kernel,tvec,imethod,0);
p_atdata(p_atdata<1e-20)=1e-20;
%use=find(p_atdata>1e-20); ntot=size(use,1);

Gp_atdata=interp1(tgrid,Gp1,tvec,imethod,0);  % for the Gaussian kernel
Gp_atdata(Gp_atdata<1e-20)=1e-20;
%Guse=find(Gp_atdata>1e-20); Gntot=size(Guse,1);

%  Find the adaptive kernel bandwidth based on the Gaussian convolution:
G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-alpha;
scale_1t = max(width/1000, min(1e10 , h_0.*lambda));

GG = exp(sum(log(Gp_atdata)) / ndata); % sum over all of the estimated density
Glambda = (Gp_atdata / GG).^-Galpha;
Gscale_1t = max(width/1000, min(4e10 , Gh_0.*Glambda));  %Gh_0 is unchanged

%%%%%%  Mike - this one loop is probably the slowest bit. Each p1new(j) does not
%%%%%%  depend on others, so if timesneeded and kernelinterp were not
%%%%%%  re-used then this could be parallel looped.

parfor j = 1:ninterp
    tnow=tgrid(j);
%    timesneeded=center+(tnow-tvec); 
%    timesneeded=center+(tnow-tvec).*(width./scale_1t).^1.0; 
    timesneeded=m_1+(tnow-tvec).*(width./scale_1t).^1.0; 
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0) ;
%    p1new(j)=(1/mass)*sum( npart.*kernelinterp );
    p1new(j)=(1/mass)*sum( npart.*(width./scale_1t).*kernelinterp );
end

if iter<4   % Gaussian always converges fast
 parfor j = 1:ninterp
    tnow=tgrid(j);
%    timesneeded=center+(tnow-tvec).*1./Gscale_1t.^1.0; 
    timesneeded=center+(tnow-tvec).*1./Gscale_1t.^1.0; 
    Gkernelinterp=interp1(tgrid,Gkernel,timesneeded,imethod,0) ;
    Gp1new(j)=(1/mass)*sum( npart.*(1./Gscale_1t.^1.0).*Gkernelinterp );
 end
end
%apple=find(p1new>1e-10); nuse=length(apple);
%rmsenow=sqrt(  dt*sum( (p1new(apple)-p1(apple)).^2)  )

%blah=((p1new-p1).^2)';
rmsenow=sqrt(  trapz(tgrid, (p1new-p1).^2) )
p1=p1new;
kernel=p1new;
intdens=trapz(tgrid,kernel)
%intdens=0.5*sum(dt.*(kernel(1:end-1)+kernel(2:end))')
kernel=kernel./intdens;

Gp1=Gp1new';

%if(rmsenow>rmselast)
%    h_0=0.8*h_0
%end
%rmselast=rmsenow;

if(rmsenow>rmselast)
    h_0=max([1*width*ndata^(-1/3) 0.9*h_0]);
    rmselast=1;  % Gives at least one iteration on new h_0
else
    rmselast=rmsenow;
end
h_over_width=h_0/width

%  Make some plots of the progress ...
figure(1)
semilogy(tvec, npart, 'k+','LineWidth',.1)
hold on

semilogy(tgrid, p1, 'ko')
semilogy(tgrid, Gp1, 'r','LineWidth',1)
semilogy(tgrid, p_Gauss, 'g')
%plot(tgrid,density(mu,sig,tgrid),'b--')
%axis([-10000 1e6 1e-11 1e-3])  % For the BTC problem
axis([min(tgrid) max(tgrid) 0.0000001/ndata 1.5*max(p1)]);  % For a general problem
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel','Init. Gauss','Reality','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(2)
semilogx(tvec, 1.5*max(p1)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on

plot(tgrid, p1, 'ko')
plot(tgrid, Gp1, 'r','LineWidth',2)
plot(tgrid, p_Gauss, 'g')
%plot(tgrid,density(mu,sig,tgrid),'b--')
%axis([1 1e6 0 1e-3])   % For the BTC problem
axis([min(tgrid) max(tgrid) 0 1.5*max(p1)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel','Init. Gauss','Reality','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(3)
plot(tvec, 1.5*max(p1)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on

plot(tgrid, p1, 'ko')
plot(tgrid, Gp1, 'r','LineWidth',2)
plot(tgrid, p_Gauss, 'g')
%plot(tgrid,density(mu,sig,tgrid),'b--')
%axis([1 1e6 0 1e-3])   % For the BTC problem
axis([min(tgrid) max(tgrid) 0 1.5*max(p1)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel','Init. Gauss','Reality','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

end   % kernel iteration loop .........................................
if iter==itermax
    fail=fail+1;
%    continue
end

end   % ensemble loop


%*************
% Plot results
%*************

figure(3)
plot(tvec, npart, '.')
hold on
plot(tgrid, p_0t, 'k')
plot(tgrid, p_Gauss, 'r')
%plot(tgrid,reald,'b--')
title('Pre-iterated, Gaussian estimates');
axis([min(tgrid) max(tgrid) 0 1.5*max(p1)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
axis([min(tgrid) max(tgrid) 0 1.5*max(p1)]);
legend('Data','constant bandwidth','variable bandwidth','Histogram')
xlabel('Data Values')
ylabel('PDFs')
hold off

figure(4)
plot(tvec, 1e-3*ones(size(tvec)), '+')
hold on
plot(tgrid, p1, 'k')
plot(tgrid, Gp1, 'r')
%plot(tgrid,reald,'b--')
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
axis([min(tvec) max(tvec) 0 1.5*max(p_0t)]);
legend('Data','data kernel','Gauss kernel','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
%plot(tvec,kernel,'-r')
hold off

% Show the kernel for a few of the points
times2plot=[2 floor(ndata/2) ndata-1];
for j = 1:length(times2plot)
    nposition=times2plot(j);
    tnow=tvec(nposition)
    figure(6);
    semilogy(tnow+(tgrid-center).*(scale_1t(nposition)/(width)),((width)./scale_1t(nposition)).*kernel,'b');
    hold on;
    plot(tnow+(tgrid-center).*Gscale_1t(nposition),((1/Gscale_1t(nposition))).*Gkernel,'r');
    plot(tnow,0.01,'d');
end
plot(tvec(npart>0.5), 0.01*npart(npart>0), '.')
title('Kernels at three data points (early, middle, late)')
axis([min(tgrid) max(tgrid) 1e-10 10]);
hold off

figure(5); hold off;
% Output a subset of the smoothed BTC:
%smoothbtc = p_1t(51:460);
%save('smoothbtc.txt','-ascii','smoothbtc')

%exit
%  Might need to put some multi-line matlab functions here to generate the
%  densities of discontinuous functions.

 function y = expdens(mu,sig,x)
     y=zeros(size(x));
     y(x>=mu-1/sig)=(1/sig)*exp(-(x(x>=mu-1/sig)-mu+1/sig)/sig);
 end
  