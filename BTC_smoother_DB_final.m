clear all
disp('Beginning smoother...');

% Import times and densities from PT simulation
% fid   = load('well.extraction.txt');
% Define some variables
alpha   = 0.5;     % Adaptive factor exponent
Galpha  = 0.5;     %    " for Gauss kernel
ninterp = 10000;    % Number of linearly spaced points to get density estimate
imethod = 'linear' % Matlab's interpolation scheme (linear or spline)
nens    = 10;      % Number in ensemble (choose 1 for BTC data)
itermax = 300;     % Max number of iterations
tol     = 1e-9;    % Inter-iteration density RMSE tolerance 
ndata   = 1000;    % Number of "fake" data to generate
mu=100; sig=10;       % mean and scale of data to be generated

% Placeholders
dkl2=0; Gl2=0; G1passl2=0; Gh0l2=0;
dkl2par=0; Gl2par=0; G1passl2par=0; Gh0l2par=0;
PMISE=0; DBMISE=0; PFGMISE=0; DBh2MISE=0; DBWMISE=0; PFGh1MISE=0;

figures=1;  % change to 1 to get figures drawn during iteration

for l=1:nens    % Ensemble loop

%  Generate (or read) some data after choosing type 
%rv='IG';
%rv='Gauss';
%rv='Cauchy';
rv='Stable';
%rv='Exponential';
%rv='Levy';
%rv='readfile'   % Choose this to read a file of data, like BTC data.
if(strcmpi(rv,'readfile'))
      ranvec=load('BTC_50000');
      ranvec=ranvec(ranvec(:,1)>1e-20,1);
      ndata=length(ranvec);
      density =@(mu,sig,tvec) zeros(size(tvec));
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
end
%  Important: include all real points and interpolated (grid) points when making calculations
tgrid=union(tgrid, tvec');
ninterp=length(tgrid);
dt=diff(tgrid);
realdens=density(mu,sig,tgrid);
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

% Get some statistics using raw data not binned (i.e., raw particle arrival times)
apple=find(npart>0);
mass=sum(npart);
center=sum(npart(apple).*tvec(apple))/mass
datavariance=sum((tvec(apple)-center).^2)/mass
iqr=tvec(apple(floor(0.75*ndata))) - tvec(apple(floor(0.25*ndata)))  %interquartile range
width=iqr/1.35;
width=iqr/1.1;  % Exponential
width=iqr/2.0;  % Cauchy
width=iqr/1.5  % Sort of average?

% Choose a method to estimate initial global bandwidth h_0
h_0FT = FT_h(tvec)              % 1) Use Fourier transform routine
%h_0 = 1*sig*1.06*ndata^-0.2     % 2) Taylor series approx.
%[h0an,h0num] = Gauss_h_0 (sqrt(datavariance),ndata);  3) assume Gaussian data
%h_0 = h0num     % 3a) Exact for Gaussian kernel and a realiz. of Gauss data 
%h_0 = h0an      % 3b) Exact for Gauss kernel AND exactly Gauss data
h_0_plug=h_plug_in(tvec)  % 4) Plug-in method of Engel et al. 
h_0=h_0_plug;
h_0=h_0FT;
%h_0=10         % Manually put in whatever h_0 you want 
Gh_0=h_0;       % Use a separate h_0 for the Gaussian kernel

%%%%%%%%%%%%%%%%%%%%%%%   Now make some density interpolations %%%%%%%%%%%%%%%%%%%%%

%  First try, using Gaussian kernel and global bandwidth ...
for j = 1:ninterp
    tnow=tgrid(j);
    p_0t(j) =(1/mass)*sum( npart./sqrt(2*pi*(h_0)^2) .* exp( -(tvec-tnow).^2 / (2*(h_0)^2)) );
end

% The adaptive scheme now uses to  density estimate to get h at each data point 
% Need the value of the interpolated density at each data point - might 
% need interpolate from p_0t if grid not equal to data points
ECDF = (1:length(tvec))'/length(tvec);

%  Find the adaptive kernel bandwidth at each data point based on the Gaussian convolution:

p_atdata=interp1(tgrid,p_0t,tvec,imethod,0);
p_atdata(p_atdata<1e-50)=1e-50;    % If spline is used, can be negative density values
G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-Galpha;
h_1t = h_0.*lambda;  % may need to put real upper bound

% Make Pedretti's UAB h_2 using pre-sorted arrival times:
h_2  = h_0*(1-ECDF) + ECDF.*h_1t;
% Standard Gaussian kernel to start, using same notation as when I switch to kernel
% based on p(t):

%Get Pedretti's results  
for j=1:ninterp    % size of my interpolation grid
    tnow=tgrid(j);
    p_1t(j) = (1/mass)*sum( (npart./sqrt(2*pi*h_1t.^2)) .* exp(-(tvec-tnow).^2 ./ (2*h_1t.^2)) );  
    p2(j)   = (1/mass)*sum( (npart./sqrt(2*pi*h_2.^2) ) .* exp(-(tvec-tnow).^2 ./ (2*h_2.^2)) );  
end

if figures>0
figure(1)
semilogy(tvec, npart, 'k+','LineWidth',.1)
hold on
%semilogy(tgrid,p_0t)
semilogy(tgrid, p_1t, 'g','LineWidth',2)
semilogy(tgrid, p2, 'r','LineWidth',2)
%semilogy(tgrid, p_Gauss, 'g')
plot(tgrid,realdens,'b--','Linewidth',2)
axis([min(tgrid) max(tgrid) 0.0000001/ndata 1.5*max(p_0t)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','h_1','h_2','Reality','Histogram')
title('Gauss kernel, UAB')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(2)
loglog(tvec, 1.5*max(p_1t)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on
plot(tgrid, p_1t, 'g','Linewidth',2)
plot(tgrid, p2, 'r','LineWidth',2)
%plot(tgrid, p_Gauss, 'g')
plot(tgrid,realdens,'b--','Linewidth',2)
axis([min(tvec) max(tvec) 1e-10 1.5*max(p_0t)]);
%axis([1 1e6 1e-11 1e-2]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','h_1','h_2','Reality','Histogram')
title('Gauss kernel, UAB')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(3)
semilogx(tvec, 1.5*max(p1)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on

plot(tgrid, p_1t, 'k','LineWidth',2)
plot(tgrid, p2, 'r','LineWidth',2)
plot(tgrid,realdens,'b--','LineWidth',2)
%axis([1 1e6 0 1e-3])   % For the BTC problem
axis([min(tgrid) max(tgrid) 0 1.5*max(p_0t)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','h_1','h_2','Reality','Histogram')
title('Gauss kernel, UAB')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off
end  % draw figures


% Store the mean MISE for known densities
MISE=trapz(tgrid,(p2-realdens).^2);
PFGMISE = ((l-1)*PFGMISE+MISE)/l   % Pedretti's mean MISE w/h2
h1MISE=trapz(tgrid,(p_1t'-realdens).^2);
PFGh1MISE = ((l-1)*PFGh1MISE+h1MISE)/l   % Pedretti's mean MISE w/h1
h_0P=h_0;

%%%%%%%%%%%%%%%%%%%%%%%   Now iterate the kernel %%%%%%%%%%%%%%%%%%%%%
factor=1.05;
%h_0=h_0FT;
h_0=h_0_plug
rmselow=1000;
%h_0=h_0/5;
% The adaptive scheme now uses to  density estimate to get h at each data point 
% Need the value of the interpolated density at each data point - might 
% need interpolate from p_0t if grid not equal to data points

for j = 1:ninterp
    tnow=tgrid(j);
    p_0t(j) =(1/mass)*sum( npart./sqrt(2*pi*(h_0)^2) .* exp( -(tvec-tnow).^2 / (2*(h_0)^2)) );
end

p_atdata=interp1(tgrid,p_0t,tvec,imethod,0);
p_atdata(p_atdata<1e-20)=1e-20;    % If spline is used, can be negative density values

%  Find the adaptive kernel bandwidth at each data point based on the Gaussian convolution:
G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-Galpha;
h_1t = max(width/1000, min(1000*width, h_0.*lambda) );  %may need to put real upper bound
h_2  = h_0*(1-ECDF) + ECDF.*h_1t;
% Standard Gaussian kernel to start, using same notation as when I switch to kernel
% based on p(t):
center=0;
kernel=(1/sqrt(2*pi*1^2))*exp((tgrid-center).^2/(-2*1^2));  

parfor j=1:ninterp    % size of my interpolation grid
    tnow=tgrid(j);
    timesneeded=center+(tnow-tvec)./h_1t.^1.0;     %scaled times from current grid point to each data point
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0);
    p_1t(j)=(1/mass)*sum( npart.*(1./h_1t.^1.0).*kernelinterp );
end

p_Gauss=p_1t;  % without a better estimate of h, this is the best I can do.

%DB: Assume that p_1t is a good first estimate of the BTC, using Gaussian kernels.
% Now use the BTC itself as the kernel.  

iqrvec=quantile(tvec,[0.25 0.75]);
iqr=iqrvec(2)-iqrvec(1);
width=iqr/1.1;   % Exponential
width=iqr/2.0;   % Cauchy
width=iqr/1.5   % Sort of average?

kernel=p_1t;  % use the first estimate of BTC as the kernel
intdens=trapz(tgrid,kernel);

Gkernel=(1/sqrt(2*pi))*exp((tgrid-center).^2/(-2));

parfor j = 1:ninterp
    tnow=tgrid(j);
    timesneeded=center+(tnow-tvec).*width./h_0; 
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0) ;
    p1(j)=(1/mass)*sum( npart.*(width./h_0^1.0).*kernelinterp );
end

h2=p1;
%   Iterate kernel shape and size ........................................

rmselast=1e6; iter=0; fail=0;
kernel=max(0,p1);
kernel=kernel./trapz(tgrid,kernel);   % Normalize kernel to unit area
h2kernel=kernel;
while rmselast>tol&iter<itermax
%  recalculate the kernel scaling factor

m_0=trapz(tgrid,kernel);   %check for unit mass of density
cum=cumsum(0.5*(dt.*(kernel(1:end-1)+kernel(2:end))'));
m_1=trapz(tgrid,(tgrid.*kernel')) ;

q3=min(find(cum>0.75));
q1=max(find(cum<0.25));
iqr=(tgrid(q3)-tgrid(q1));

width=iqr/1.35;   % Gaussian
width=iqr/1.1;    % Exponential
width=iqr/2.0;    % Cauchy
width=iqr/1.5    % Sort of average?
iter=iter+1

%  Find the adaptive kernel bandwidth based on the Gaussian convolution:

p_atdata=interp1(tgrid,kernel,tvec,imethod,0);
p_atdata(p_atdata<1e-50)=1e-50;
G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-alpha;
%scale_1t = max(width/1000, min(1e10 , h_0.*lambda));
scale_1t=h_0.*lambda;
h_2  = h_0*(1-ECDF) + ECDF.*scale_1t;   % Try Pedretti's

%%%%%%  Each p1new(j) does not depend on others, so if timesneeded and 
%%%%%%  kernelinterp were not re-used then these could be parallel looped.

parfor j = 1:ninterp
    tnow=tgrid(j);
    timesneeded=m_1+(tnow-tvec).*(width./scale_1t); 
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0) ;
    p1new(j)=(1/mass)*sum( npart.*(width./scale_1t).*kernelinterp );
end

 parfor j = 1:ninterp
    tnow=tgrid(j);
    timesneeded=m_1+(tnow-tvec).*(width./h_2); 
    kernelinterp=interp1(tgrid,kernel,timesneeded,imethod,0) ;
    h2new(j)=(1/mass)*sum( npart.*(width./h_2).*kernelinterp );
 end
rmsenow=sqrt(  trapz(tgrid, (p1new-p1).^2) )
p1=p1new;
h2=h2new';
kernel=p1new;

intdens=trapz(tgrid,kernel)
kernel=kernel./intdens;

if(rmsenow<rmselow) ; rmselow=rmsenow; end
if(rmselow<1e-4);     factor=1.1   ; end
if(rmselow<1e-5);     factor=1.05    ; end
if(rmselow<1e-6);     factor=1.02    ; end

if(rmsenow>rmselast || intdens <0.98)
    h_0=factor*h_0
    rmselast=1;  % Gives at least one iteration on new h_0
else
    rmselast=rmsenow;
end

%  Make some plots of the progress ...............
if figures>0
figure(11)
semilogy(tvec, npart, 'k+','LineWidth',.1)
hold on

semilogy(tgrid, p1, 'ko')
semilogy(tgrid, h2, 'r','LineWidth',1)
semilogy(tgrid, p_Gauss, 'g')
plot(tgrid,realdens,'b--')
%axis([-10000 1e6 1e-11 1e-3])  % For the BTC problem
axis([min(tgrid) max(tgrid) 0.0000001/ndata 1.5*max(p1)]);  % For a general problem
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel+h2','Init. Gauss','Reality','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(13)
semilogx(tvec, 1.5*max(p1)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on

plot(tgrid, p1, 'ko')
plot(tgrid, h2, 'r','LineWidth',2)
plot(tgrid, p_Gauss, 'g')
plot(tgrid,realdens,'b--')
%axis([1 1e6 0 1e-3])   % For the BTC problem
axis([min(tgrid) max(tgrid) 0 1.5*max(p1)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel','Init. Gauss','Reality','Histogram')
title('Iterated kernels')
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off

figure(12)
loglog(tvec, 1.5*max(p1)*ones(size(tvec)), 'k+','LineWidth',.1)
hold on

plot(tgrid, p1, 'ko')
plot(tgrid, h2, 'r','LineWidth',2)
%plot(tgrid, p_Gauss, 'g')
%plot(tgrid,realdens,'b--')
axis([1 1e6 1e-11 1e-2])   % For the BTC problem
%axis([min(tvec) max(tvec) 0 1.5*max(p1)]);
histogram(tvec(npart>0),'Normalization','pdf','FaceAlpha',0.05,'NumBins',50)
legend('Data','Data kernel','Gauss kernel','Init. Gauss','Reality','Histogram')
title(['Iterated kernels, h_0=',num2str(h_0)])
xlabel('Data Values')
ylabel('PDFs')
drawnow
hold off
end % if figures

end   % kernel iteration loop .........................................
if iter==itermax
    fail=fail+1;
%    continue
end
MISE=trapz(tgrid,(p1'-realdens).^2);
DBMISE = ((l-1)*DBMISE+MISE)/l   % Benson's mean MISE
WMISE=trapz(tgrid,(1./p1').*(p1'-realdens).^2);
DBWMISE = ((l-1)*DBWMISE+WMISE)/l   % Benson's mean weighted MISE
h2MISE=trapz(tgrid,(h2'-realdens).^2);
DBh2MISE = ((l-1)*DBh2MISE+h2MISE)/l   % Benson's mean MISE using h2
h_0DB=h_0
pause(2)
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
plot(tgrid, h2, 'r')
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
%    plot(tnow+(tgrid-center).*Gscale_1t(nposition),((1/Gscale_1t(nposition))).*Gkernel,'r');
    plot(tnow,0.01,'d');
end
plot(tvec(npart>0.5), 0.01*npart(npart>0), '.')
title('Kernels at three data points (early, middle, late)')
axis([min(tgrid) max(tgrid) 1e-10 10]);
hold off

figure(5); hold off;

%  Might need to put some multi-line matlab functions here to generate the
%  densities of discontinuous functions.

 function y = expdens(mu,sig,x)
     y=zeros(size(x));
     y(x>=mu-1/sig)=(1/sig)*exp(-(x(x>=mu-1/sig)-mu+1/sig)/sig);
 end
  