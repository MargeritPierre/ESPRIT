%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SNR)

clc
clear all
%close all


Ns = 50 ; % Number of Samples
F = [-.023 .1 -.2 .06 .3]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = logspace(-1,4,10) ;
nMCMC = 200 ;
profiler = false ;

t = 0:1:Ns-1 ; % time
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

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
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1))*U ;
        Signal = sum(Amp.'*sum(exp(1i*F.'*t),1),1) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        [~,ind] = sort(real(out.K)) ;
        K(s,:,m) = out.K(ind) ;
        dK(s,:,m) = out.dK(ind) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nSNR,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

dK = reshape(dK,nSNR,[]) ;
stdK = mean(var(K,0,3),2) ;
meandK = mean(abs(dK),2) ;
mindK = min(abs(dK),[],2) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],2) ; ... meandK + std(dK,0,3) ;

[~,ind] = sort(real(F)) ;
SE = abs(K-repmat(sort(F(ind)),[nSNR 1 nMCMC])).^2 ;
SE = reshape(SE,nSNR,[]) ;
tMSE = (mean(SE,2)) ;
minSE = min(SE,[],2) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],2) ; ... meandK + std(dK,0,3) ;


clf ;
plot(SNR,mindK,':k') ;
plot(SNR,meandK,'.-k','markersize',20) ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
plot(SNR,maxSE,':r') ;
plot(SNR,minSE,':r') ;
plot(SNR,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SIGNAL SIZE)

clc
clear all
%close all


Ns = round(logspace(log10(5),log10(200),20)) ; % Number of Samples
F = [-.23-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 1e1 ;
nMCMC = 100 ;
profiler = false ;

nF = length(F) ; % number of tones
nNs = length(Ns) ; % number of signal sizes

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

K = zeros(nNs,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for n = 1:nNs
    t = 0:1:Ns(n)-1 ; % time 
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1))*U ;
        Signal = sum(Amp.'*sum(exp(1i*F.'*t),1),1) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(n,:,m) = out.K ;
        dK(n,:,m) = out.dK ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(n-1)*nMCMC)/nMCMC/nNs,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

meanK = mean(K,3) ;
stdK = var(K,0,3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;
SE = abs(K-repmat(F.',[nNs 1 nMCMC])).^2 ;
tMSE = (mean(SE,3)) ;
minSE = min(SE,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],3) ; ... meandK + std(dK,0,3) ;


clf ;
plot(Ns,mindK,':k') ;
plot(Ns,meandK,'.-k','markersize',20) ;
plot(Ns,maxdK,':k') ;
plot(Ns,tMSE,'.r','markersize',35) ;
plot(Ns,maxSE,':r') ;
plot(Ns,minSE,':r') ;
plot(Ns,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE POLE FREQUENCY)

clc
clear all
%close all


Ns = 100 ; % Number of Samples
F = logspace(log10(0.000001),log10(Ns/2.1),10)*pi/Ns ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 1e4 ;
nMCMC = 100 ;
profiler = false ;

nF = length(F) ; % number of tones
t = 0:1:Ns-1 ; % time 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , 1 ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
          }' ;

K = zeros(nF,1,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for f = 1:nF
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1))*U ;
        Signal = sum(U.'*sum(exp(1i*F(f).'*t),1),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(f,:,m) = out.K ;
        dK(f,:,m) = out.dK ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(f-1)*nMCMC)/nMCMC/nF,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

stdK = var(K,0,3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;
tMSE = (mean(abs(K-repmat(F.',[1 1 nMCMC])).^2,3)) ;


clf ;
plot(F,mindK,':k') ;
plot(F,meandK,'.-k','markersize',20) ;
plot(F,maxdK,':k') ;
plot(F,tMSE,'.r','markersize',35) ;
plot(F,stdK,'ob','markersize',20) ;
set(gca,'xscale','log') ;
set(gca,'yscale','log') ;



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SNR, Multiple Snapshot)

clc
clear all
%close all


Ns = 100 ; % Number of Sensors
Nsnap = 100 ; % Number of Snapshots
F = [-.023-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = logspace(-1,4,10) ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels
Amp = ((rand(Nsnap,1)*2-1)+1i*(rand(Nsnap,1)*2-1))*U ;
Signal = Amp*sum(exp(1i*F.'*t),1) ; 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [Ns] ; ... 
          }' ;

K = zeros(nSNR,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
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

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

meanK = mean(K,3) ;
stdK = var(K,0,3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;
tMSE = (mean(abs(K-repmat(F.',[nSNR 1 nMCMC])).^2,3)) ;


clf ;
plot(SNR,mindK,':k') ;
plot(SNR,meandK,'.-k','markersize',20) ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
plot(SNR,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE NUMBER OF SNAPSHOTS)

clc
clear all
%close all


Ns = 20 ; % Number of Sensors
Nsnap = unique(round(logspace(log10(1),log10(500),10))) ; % Number of Snapshots
F = [-.23-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 1e2 ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nNsnap = length(Nsnap) ; % number of number of snapshots

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [Ns] ; ... 
          }' ;

K = zeros(nNsnap,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nNsnap
    for m = 1:nMCMC
        Amp = ((rand(Nsnap(s),1)*2-1)+1i*(rand(Nsnap(s),1)*2-1))*U ;
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(s,:,m) = out.K ;
        dK(s,:,m) = out.dK ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nNsnap,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

meanK = mean(K,3) ;
stdK = var(K,0,3) ;

meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;

SE = abs(K-repmat(F.',[nNsnap 1 nMCMC])).^2 ;
tMSE = (mean(SE,3)) ;
minSE = min(SE,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],3) ; ... meandK + std(dK,0,3) ;

clf ;
plot(Nsnap,mindK,':k') ;
plot(Nsnap,meandK,'.-k','markersize',20) ;
plot(Nsnap,maxdK,':k') ;
plot(Nsnap,tMSE,'.r','markersize',35) ;
plot(Nsnap,maxSE,':r') ;
plot(Nsnap,minSE,':r') ;
plot(Nsnap,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;








%% 1D UNCERTAINTY ESTIMATION (VARAIATION OF THE SPATIAL SMOOTHING)

clc
clear all
%close all


Ns = 100 ; % Number of Samples
F = [-.023-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
%K0 = unique(round(logspace(log10(3),log10(Ns),10))) ;
K0 = unique(round(linspace(3,Ns,20))) ;
SNR = 1e2 ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nK0 = length(K0) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
          }' ;

K = zeros(nK0,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for k = 1:nK0
    ARGS = [EspArgs{:},{'K0'},{K0(k)}] ;
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1))*U ;
        Signal = sum(Amp.'*sum(exp(1i*F.'*t),1),1) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,ARGS{:}) ;
        K(k,:,m) = out.K ;
        dK(k,:,m) = out.dK ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(k-1)*nMCMC)/nMCMC/nK0,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

%[K,ind] = sort(K,2) ;
%dK = dK(ind) ;

meanK = mean(K,3) ;
stdK = var(K,0,3) ;
meandK = mean(abs(dK),3) ;
mindK = min(abs(dK),[],3) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],3) ; ... meandK + std(dK,0,3) ;
SE = abs(K-repmat(F.',[nK0 1 nMCMC])).^2 ;
tMSE = (mean(SE,3)) ;
minSE = min(SE,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],3) ; ... meandK + std(dK,0,3) ;


clf ;
plot(K0,mindK,':k') ;
plot(K0,meandK,'.-k','markersize',20) ;
plot(K0,maxdK,':k') ;
plot(K0,tMSE,'.r','markersize',35) ;
plot(K0,maxSE,':r') ;
plot(K0,minSE,':r') ;
plot(K0,stdK,'ob','markersize',20) ;
%set(gca,'xscale','log') ;
set(gca,'yscale','log') ;




