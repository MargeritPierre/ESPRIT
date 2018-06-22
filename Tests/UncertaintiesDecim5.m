%% PROBABILITY OF DETECTION AS A FUNCTION OF SNR

clc
clear all
%close all

%clf ;

Nsens = 100 ; % Number of Sensors
Nsnap = 1 ; % Number of Snapshots
F = [-0.05 0.04] ; 0+[-1 1]*0.05 ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = logspace(-2,6,10) ;
nMCMC = 100 ;

t = 0:1:Nsens-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , 1:10 ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [6] ; ...
           'COMPUTE_U' , false ; ...
           'COMPUTE_dK' , false ; ...
           'CRITERION' , 'MDL' ; ...
           'CRIT_THRS' , 1 ; ...
          }' ;

Ndetect = zeros(nSNR,1,nMCMC) ;
wtbr = waitbar(0) ;
ti = tic ;
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(Nsnap,1)*2-1)+1i*(rand(Nsnap,1)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        Ndetect(s,:,m) = length(out.K) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nSNR,wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;

Detected = Ndetect==nF ;
probDetect = mean(Detected,3) ;

plot(SNR,probDetect,'k') ;
set(gca,'xscale','log') ;
%set(gca,'yscale','log') ;


%% PROBABILITY OF DETECTION AS A FUNCTION OF DECIMATION FACTOR

clc
clear all
%close all

clf ;

Nsens = 100 ; % Number of Sensors
Nsnap = 1 ; % Number of Snapshots
F = [-0.05 0.04] ; 0+[-1 1]*0.05 ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 10^.8 ;
nMCMC = 100 ;
decim = 1:10 ;

t = 0:1:Nsens-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           ... 'R0' , 1:10 ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [6] ; ...
           'COMPUTE_U' , false ; ...
           'COMPUTE_dK' , false ; ...
           'CRITERION' , 'ESTER' ; ...
           'CRIT_THRS' , 1 ; ...
          }' ;

Ndetect = zeros(length(decim),1,nMCMC) ;
wtbr = waitbar(0) ;
ti = tic ;
for d = 1:length(decim)
    for m = 1:nMCMC
        Amp = ((rand(Nsnap,1)*2-1)+1i*(rand(Nsnap,1)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)^2/norm(noise)^2/SNR ;
        Args = [EspArgs,[{'DECIM'};{[1 decim(d)]}]] ;
        out = ESPRIT(Signal+noise,Args{:}) ;
        Ndetect(d,:,m) = length(out.K) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(d-1)*nMCMC)/nMCMC/length(decim),wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;

Detected = Ndetect==nF ;
probDetect = mean(Detected,3) ;

plot(decim,probDetect,'k') ;
%set(gca,'xscale','log') ;
%set(gca,'yscale','log') ;


%% VARIANCE AS A FUNCTION OF DECIMATION FACTOR

clc
clear all
%close all

clf ;

Nsens = 50 ; % Number of Sensors
Nsnap = 1 ; % Number of Snapshots
F = [-0.05 0.04]/2 ; 0+[-1 1]*0.05 ; % Tones (normalized freq.)
U = [1].' ; % Amplitudes
SNR = 10^.7 ;
nMCMC = 100 ;
decim = 1:10 ;

t = 0:1:Nsens-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
           ...'COMPUTE_U' , false ; ...
           ...'COMPUTE_dK' , false ; ...
           'CRITERION' , 'ESTER' ; ...
           'CRIT_THRS' , 1 ; ...
          }' ;

K = zeros(length(decim),nF,nMCMC) ;
dK = zeros(length(decim),nF,nMCMC) ;
wtbr = waitbar(0) ;
ti = tic ;
for d = 1:length(decim)
    for m = 1:nMCMC
        Amp = ((rand(Nsnap,1)*2-1)+1i*(rand(Nsnap,1)*2-1)) ;
        Amp = Amp./abs(Amp).*U ;
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        Args = [EspArgs,[{'DECIM'};{[1 decim(d)]}]] ;
        out = ESPRIT(Signal+noise,Args{:}) ;
        [kkk,ind] = sort(real(out.K)) ;
        K(d,:,m) = kkk ;
        dK(d,:,m) = out.dK(ind) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(d-1)*nMCMC)/nMCMC/length(decim),wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;

varK = var(K,0,3) ;

plot(decim,varK,'or') ;
plot(decim,mean(abs(dK),3),'k') ; 
%set(gca,'xscale','log') ;
set(gca,'yscale','log') ;
 
 
 
%% VARIANCE AS A FUNCTION OF THE SIGNAL MESH STEP

clc
clear all
%close all

clf ;

L = 100;
Nsens = unique(round(linspace(10,100,10))) ; % Number of Sensors
Nsnap = 1 ; % Number of Snapshots
F = [-0.05 0.04] ; 0+[-1 1]*0.05 ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 10^2 ;
nMCMC = 100 ;
nF = length(F) ; % number of tones

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           ...'K0' , [] ; ...
           'COMPUTE_U' , false ; ...
           'COMPUTE_dK' , false ; ...
           'CRITERION' , 'MDL' ; ...
           'CRIT_THRS' , 1 ; ...
          }' ;

K = zeros(length(Nsens),nF,nMCMC) ;
wtbr = waitbar(0) ;
ti = tic ;
for d = 1:length(Nsens)
    t = linspace(0,L,Nsens(d)) ; % time 
    for m = 1:nMCMC
        Amp = ((rand(Nsnap,1)*2-1)+1i*(rand(Nsnap,1)*2-1)) ;
        Amp = repmat(U,[Nsnap 1]) ; Amp./abs(Amp)*U ;
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)^2/norm(noise)^2/SNR ;
        Args = [EspArgs,[{'K0'};{Nsens(d)}]] ;
        out = ESPRIT(Signal+noise,Args{:}) ;
        K(d,:,m) = sort(real(out.K)) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(d-1)*nMCMC)/nMCMC/length(Nsens),wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;

varK = var(K,0,3) ;

plot(Nsens,varK,'k') ;
%set(gca,'xscale','log') ;
 set(gca,'yscale','log') ;







