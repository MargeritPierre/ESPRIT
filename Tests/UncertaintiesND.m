%% ND UNCERTAINTY ESTIMATION (VARIATION OF THE SNR)

clc
clear all
%close all


Lk = [10 10 10] ;
F = [.03 ; .04 ; -.02]*pi ; % .1 -.2 .06 .3 -.04 .01]*pi ; % Tones (normalized freq.)];%
U = [1] ; % Amplitudes
SNR = logspace(-1,6,10) ;
nMCMC = 100 ;
profiler = false ;

xx = 0:1:Lk(1)-1 ; yy = 0:1:Lk(1)-1 ; zz = 0:1:Lk(1)-1 ;
[XX,YY,ZZ] = meshgrid(xx,yy,zz) ;
nF = size(F,2) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , [1 2 3] ; ...
           'R0' , nF ; ...
           'DECIM' , ones(size(Lk)) ; ...
           'FUNC' , repmat('exp',[length(Lk) 1]) ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
          }' ;

K = zeros(nSNR,length(Lk),nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        Signal = Amp*exp(1i*(XX*F(2,:)+YY*F(1,:)+ZZ*F(3,:))) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal(:))/norm(noise(:))/SNR(s) ;
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
if profiler ; profile viewer ; end
drawnow ;

dK = reshape(dK,nSNR,length(Lk),[]) ;
stdK = mean(var(K,0,3),3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;

SE = abs(K-repmat(F.',[nSNR 1 nMCMC])).^2 ;
SE = reshape(SE,nSNR,length(Lk),[]) ;
tMSE = (mean(SE,3)) ;
minSE = min(SE,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],3) ; ... meandK + std(dK,0,3) ;


clf ;
plot(SNR,mindK,':k') ;
plot(SNR,meandK,'.-k','markersize',20) ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
plot(SNR,maxSE,':r') ;
plot(SNR,minSE,':r') ;
plot(SNR,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;



