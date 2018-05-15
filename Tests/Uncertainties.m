%% 1D UNCERTAINTY ESTIMATION

clc
clear all
%close all


Ns = 100 ; % Number of Samples
F = [-.2-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = logspace(5,10,20) ;
nMCMC = 100 ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels
Signal = sum(U.'*sum(exp(1i*F.'*t),1),1) ; 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
          }' ;

K = zeros(nSNR,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
for s = 1:nSNR
    for m = 1:nMCMC
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In)*max(abs(Signal(:)))/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(s,:,m) = out.K ;
        dK(s,:,m) = out.dK ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nSNR,wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

meanK = mean(K,3) ;
stdK = var(K,0,3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ;
maxdK = max(abs(dK),[],3) ;
tMSE = sqrt(mean(abs(K-repmat(F.',[nSNR 1 nMCMC])).^2,3)) ;


clf ;
plot(SNR,mindK,':k') ;
plot(SNR,meandK,'-k') ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
set(gca,'xscale','log','yscale','log') ;








