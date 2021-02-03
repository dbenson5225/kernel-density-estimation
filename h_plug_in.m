function [h0] = h_plug_in(X)

ndata=length(X);
%mu=mean(X);
%sig=std(X);
iqrvec=quantile(X,[0.25 0.75]);  %interquartile range of data
iqr=iqrvec(2)-iqrvec(1);

RT2PI = sqrt(2.*pi);
      
% calculate global optimal h

%XIQR=X(NINT(.75*XN))-X(NINT(.25*XN)+1)
XIQR=iqr;  XN=ndata ;
ITER=5 ;
%-------  Estimate inflation constant C
H2=(.920*XIQR)/(ndata^(1./7.)) ;
H3=(.912*XIQR)/(ndata^(1./9.)) ;
S2=0; S3=0;
for I = 1:ndata-1
    for J = I+1:ndata
        D2 = ((X(I) - X(J))/H2)^2 ;
        D3 = ((X(I) - X(J))/H3)^2 ;
%        IF(D2.GT.50.AND.D3.GT.60) GOTO 20
        E2 = exp(-0.5*D2) ;
        E3 = exp(-0.5*D3) ;
        S2 = S2 + (D2^2 - 6*D2 + 3)*E2 ;
        S3 = S3 + (D3^3 - 15*D3^2 + 45*D3 - 15)*E3 ;
    end
end
RHAT2 = (2*S2)/((XN^2)*(H2^5)*RT2PI) ;
RHAT2 = RHAT2 + 3/(RT2PI*XN*(H2^5)) ;
RHAT3 = (-2*S3)/((XN^2)*(H3^7)*RT2PI) ;
RHAT3 = RHAT3 + 15./(RT2PI*XN*(H3^7)) ;
CO1 = 1.357*(RHAT2/RHAT3)^(1/7) ;
%
%-------  Compute constant of asymptotic formula
%-
   CONST=1./(2.*sqrt(pi)) ;
   A=1.132795764/RHAT3^(1./7.)*XN^(-1./2.) ;
%-
%------  Loop over iterations
%-
 for IT=1:ITER                                                                     
%-------  Estimate functional
        S=0. ;
        for I = 1:ndata
           for J = I+1:ndata
                     D2 = ((X(I) - X(J))/A)^2 ;
%                     IF(D2.GT.50) GOTO 40
                     E2 = exp(-0.5*D2) ;
                     S = S + (D2^2 - 6*D2 + 3)*E2 ;
           end
        end
        R2 = (2*S)/((XN^2)*(A^5)*RT2PI) ;
        R2 = R2 + 3./(RT2PI*XN*(A^5)) ;
%-                                                                     
%-------  Estimate bandwidth by asymptotic formula
%-
         H=(CONST/(R2*XN))^(.2) ;                            
         A=CO1*H^(5./7.) ;
 end
 
h0=H;

end
