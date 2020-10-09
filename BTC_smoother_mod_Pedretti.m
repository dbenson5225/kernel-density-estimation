clear all
disp('Beginning iterative version of Pedretti smoother...');

% Import times and densities from PT simulation
% fid   = load('well.extraction.txt');
% Define some variables
alpha   = 0.5;     % Adaptive factor exponent
Galpha  = 0.5;     %    " for Gauss kernel
ninterp = 50000;    % Number of linearly spaced points to get density estimate
imethod = 'linear' % Matlab's interpolation scheme (linear or spline)
%nens    = 10;      % Number in ensemble
itermax = 200;     % Max number of iterations
tol     = 1e-8;    % Inter-iteration density RMSE tolerance 
%ndata   = 1000;    % Number of "fake" data to generate
%mu=0; sig=1;       % mean and scale of data to be generated

% Placeholders
dkl2=0; Gl2=0; G1passl2=0; Gh0l2=0;
dkl2par=0; Gl2par=0; G1passl2par=0; Gh0l2par=0;

rv='readfile'
if(strcmpi(rv,'readfile'))
      ranvec=load('BTC_500');
      ranvec=ranvec(ranvec(:,1)>1e-20,1);
      ndata=length(ranvec);
end

fid(1:ndata,3)=1;      % This is a structure useful to importing particle arrival times
fid(1:ndata,2)=ranvec; % So I'll keep it
fid=sortrows(fid,2);   % Need to sort data so that the kernel can be interpolated later

tvec  = fid(:,2);       % X-axis values
npart = fid(:,3);       % Particle concentration (ones and zeros for observed times/padding)
ndata = size(tvec,1);   % Actual number of data 

texpand=1; tspan=max(tvec)-min(tvec);  % Need to expand range to make kernel extrapolate
tgrid=linspace(min(tvec)-(texpand*tspan),max(tvec)+1*(texpand*tspan),ninterp);  % regularly spaced points to calc. f(t).

% Manually generate interpolation grid (might be wated for BTCs)
tgrid=linspace(0,1e6,ninterp);
dt=tgrid(2)-tgrid(1)

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

apple=find(npart>0);
mass=sum(npart);
center=sum(npart(apple).*tvec(apple))/mass
datavariance=sum((tvec(apple)-center).^2)/mass
datasig=sqrt(datavariance);
iqr=1*(tvec(apple(floor(0.75*ndata))) - tvec(apple(floor(0.25*ndata))))  %interquartile range

% Use these for a first approximation:  The width of the data will be the
% "concrete" width against which the kernel is scaled.  Choose one method:

%h_0 = FT_h(tvec)               % 1) Use Fourier transform routine

[h0an,h0num] = Gauss_h_0 (sqrt(datavariance),ndata);
h_0 = h0num     % 3) Exact for Gaussian kernel and a realiz. of Gauss data 
%h_0 = h0an      % 4) Exact for Gauss kernel AND exactly Gauss data
Gh_0=h_0;       % Use a separate h_0 for the Gaussian kernel

%%%%%%%%%%%%%%%%%%%%%%%   Now make some density interpolations %%%%%%%%%%%%%%%%%%%%%

%  First try, using Gaussian kernel ...
for j = 1:ninterp
    tnow=tgrid(j);
    p_0t(j) =(1/mass)*sum( npart./sqrt(2*pi*(h_0)^2) .* exp( -(tvec-tnow).^2 / (2*(h_0)^2)) );
end

% The adaptive scheme now uses to  density estimate to get h at each data point 
% Need the value of the interpolated density at each data point - might 
% need interpolate from p_0t if grid not equal to data points
ECDF = (1:length(tvec))'/length(tvec);

%h_0=50
p2=p_0t;
for iteration=1:50
    iteration
%  Find the adaptive kernel bandwidth at each data point based on the Gaussian convolution:
%  First calculate h_0 based on density as it iterates

dk=.01;
Kgauss=1./sqrt(2*pi) .* exp(-(-5:dk:5).^2 / 2);
twodiff=(p2(1:end-2)-2*p2(2:end-1)+p2(3:end))/dt/dt;
h_0new=(sum(Kgauss.*Kgauss)/(ndata*sum(twodiff.*twodiff)*dt))^(0.2)
if(abs(h_0new - h_0)/h_0)<.001
    break   % Check convergence by h_0 change
end
h_0=min([h_0 h_0new]);

p_atdata=interp1(tgrid,p2,tvec,imethod,0);
p_atdata(p_atdata<1e-30)=1e-30;    % If spline is used, can be negative density values

G = exp(sum(log(p_atdata)) / ndata); % sum over all of the estimated density
lambda = (p_atdata / G).^-Galpha;
h_1t = max(datasig/1000, min(1000*datasig, h_0.*lambda) );  %may need to put real upper bound

% Make Pedretti's UAB h_2 using pre-sorted arrival times:
h_2  = h_0*(1-ECDF) + ECDF.*h_1t;
% Standard Gaussian kernel to start
%kernel=(1/sqrt(2*pi))*exp((tgrid-center).^2/-2);  

%Get Pedretti's first result (may need to iterate as p_atdata changes) 
parfor j=1:ninterp    % size of my interpolation grid
    tnow=tgrid(j);
    p_1t(j) = (1/mass)*sum( npart./sqrt(2*pi*h_1t.^2) .* exp(-(tvec-tnow).^2 ./ (2*h_1t.^2)) );  
    p2(j)   = (1/mass)*sum( npart./sqrt(2*pi*h_2.^2) .* exp(-(tvec-tnow).^2 ./ (2*h_2.^2)) );  
end

figure(11)
loglog(tgrid,p_0t)
axis([10 1e6 1e-10 1e-3])
hold on
plot(tgrid,p_1t)
plot(tgrid,p2)
drawnow

figure(12)
loglog(tgrid,p2)
axis([10 1e6 1e-10 1e-3])
hold on
drawnow

end

  