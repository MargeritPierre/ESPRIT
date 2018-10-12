%% INFLUENCE OF A CIRCULAR WAVE

clc
clear all

clf
kTH = 0.1*pi ;
L = [30,30] ;
x0 = [-5 -5]*3 ;

EspArgs = { ...
            'DIMS_K',[1 2] ...
            ,'R0',1:10 ...
            ,'SOLVER','eig' ...
            ,'CRITERION','ESTER' ...
            ,'CRIT_THRS', 1 ...
            ,'M/L',2/3 ...
            ,'FUNC',repmat('exp',[2 1]) ...
            ,'SHIFTS',eye(2) ...
            ,'DECIM', [1 1]*1  ...
            ,'COMPUTE_U', false  ...
            ,'COMPUTE_dK', true  ...
            ,'DEBUG', true  ...
            }.' ;

[XX,YY] = meshgrid(-floor(L(2)/2)+(1:L(2)),floor(-L(1)/2)+(1:L(1))) ;

RR = sqrt((XX-x0(2)).^2 + (YY-x0(1)).^2) ;
U = exp(1i*RR*kTH) ;

r = L ;
[u,s,v] = svd(U) ;
U=u(:,1:r)*s(1:r,1:r)*v(:,1:r)' ;

surf(XX,YY,real(U),'facecolor','interp','edgecolor','none')
axis tight
myaxisequal('xy') ;
      
OUT = ESPRIT(U,EspArgs{:}) ;
K = real(OUT.K) ;
dK = real(OUT.dK) ;
k = sqrt(sum(K.^2,1))
dk = sqrt(sum(dK.^2,1))
diffK = abs(k-kTH)







