%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SNR)

clc
clear all
%close all
%clf ;


Ns = 100 ; % Number of Samples
F = 0.01;%+0.1i ; (linspace(-1,1,1)*.04 + .03)*pi*(1-.10i) ; %rand(1,8)*0.49*pi ; %[.05 -.1 .3 -.02]*pi%] ; % Tones (normalized freq.)];%]%
U = [1] ; % Amplitudes
SNR = logspace(-2,4,10) ;
FUNC = 'exp' ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [] ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nSNR,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(1,nF)*2-1)+1i*(rand(1,nF)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        switch FUNC
            case 'exp'
                Signal = Amp*exp(1i*F.'*t) ;
            case 'cos'
                Signal = Amp*cos(F.'*t) ;
        end
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        [~,ind] = sort(real(out.K)) ;
        K(s,:,m) = out.K(ind) ;
        dK(s,:,m) = out.dK(ind).^2 ;
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

% plot(SNR,meandK,'.-','markersize',20) ;
% set(gca,'colororderindex',get(gca,'colororderindex')-1) ;
% plot(SNR,tMSE,'.','markersize',35) ;
% set(gca,'xscale','log','yscale','log') ;
% grid on


plot(SNR,mindK,':k') ;
plot(SNR,meandK,'.-k','markersize',20) ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
plot(SNR,maxSE,':r') ;
plot(SNR,minSE,':r') ;
plot(SNR,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;
grid on



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SIGNAL SIZE)

clc
clear all
%close all


Ns = round(logspace(log10(10),log10(100),20)) ; % Number of Samples
F = [.03]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
FUNC = 'cos' ;
SNR = 1e6 ;
nMCMC = 100 ;
profiler = false ;

nF = length(F) ; % number of tones
nNs = length(Ns) ; % number of signal sizes

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [1/2] ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nNs,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for n = 1:nNs
    t = 0:1:Ns(n)-1 ; % time 
    for m = 1:nMCMC
        Amp = ((rand(1,nF)*2-1)+1i*(rand(1,nF)*2-1)) ;
        Amp = 1 ; Amp./abs(Amp)*U ;
        switch FUNC
            case 'exp'
                Signal = Amp*exp(1i*F.'*t) ;
            case 'cos'
                Signal = Amp*cos(F.'*t) ;
        end
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(n,:,m) = out.K ;
        dK(n,:,m) = out.dK.^2 ;
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
F = logspace(log10(0.000001),log10(Ns/2.01),10)*pi/Ns;%*(1+.05i) ; % Tones (normalized freq.)
U = [100] ; % Amplitudes
FUNC = 'cos' ;
SNR = 1e2 ;
nMCMC = 100 ;
profiler = false ;

nF = length(F) ; % number of tones
t = 0:1:Ns-1 ; % time 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , 1 ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nF,1,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for f = 1:nF
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        switch FUNC
            case 'exp'
                Signal = Amp*exp(1i*F(f).'*t) ;
            case 'cos'
                Signal = Amp*cos(F(f).'*t) ;
        end
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(f,:,m) = out.K ;
        dK(f,:,m) = out.dK.^2 ;
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

dK = reshape(dK,nF,[]) ;
stdK = mean(var(K,0,3),2) ;
meandK = mean(abs(dK),2) ;
mindK = min(abs(dK),[],2) ; ... meandK - std(dK,0,3) ;
maxdK = max(abs(dK),[],2) ; ... meandK + std(dK,0,3) ;

[~,ind] = sort(real(F)) ;
SE = abs(K-repmat(F.',[1 1 nMCMC])).^2 ;
SE = reshape(SE,nF,[]) ;
tMSE = (mean(SE,2)) ;
minSE = min(SE,[],2) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],2) ; ... meandK + std(dK,0,3) ;


clf ;
plot(F/pi,mindK,':k') ;
plot(F/pi,meandK,'.-k','markersize',20) ;
plot(F/pi,maxdK,':k') ;
plot(F/pi,tMSE,'.r','markersize',35) ;
%plot(F/pi,maxSE,':r') ;
%plot(F/pi,minSE,':r') ;
plot(F/pi,stdK,'ob','markersize',20) ;
plot(F/pi,(F/pi).^2,'-.r','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;
grid on



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

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , 2/3 ; ... 
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nSNR,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
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


Ns = 50 ; % Number of Sensors
Nsnap = unique(round(logspace(log10(5),log10(500),10))) ; % Number of Snapshots
F = [-.0023-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = 10^2 ;
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
           'M/L' , 0 ; ... 
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nNsnap,nF,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nNsnap
    for m = 1:nMCMC
        Amp = ((rand(Nsnap(s),1)*2-1)+1i*(rand(Nsnap(s),1)*2-1)) ;
        Amp = Amp./abs(Amp).*U ;
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
%plot(Nsnap,maxSE,':r') ;
%plot(Nsnap,minSE,':r') ;
plot(Nsnap,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;








%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SPATIAL SMOOTHING)

clc
clear all
%close all


Ns = 100 ; % Number of Samples
F = 0.001-0.1i % [-.23-.1i*0]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
%K0 = unique(round(logspace(log10(3),log10(Ns),10))) ;
M_L = linspace(0.5,1,10) ; M_L = M_L(1:end-1) ; %M_L = M_L(1:end-1) ;
SNR = 1e2 ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nM_L = length(M_L) ; % number of Spatial Smoothing factors

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nM_L,nF,nMCMC) ;
dK = K ;
M_L_exp = zeros(nM_L,1) ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for k = 1:nM_L
    ARGS = [EspArgs{:},{'M/L'},{M_L(k)}] ;
    for m = 1:nMCMC
        Amp = ((rand(1)*2-1)+1i*(rand(1)*2-1)) ;
        Amp = Amp./abs(Amp)*U ;
        Signal = sum(Amp.'*sum(exp(1i*F.'*t),1),1) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,ARGS{:}) ;
        K(k,:,m) = out.K ;
        dK(k,:,m) = out.dK ;
        M_L_exp(k) = out.Mk./out.Lk ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(k-1)*nMCMC)/nMCMC/nM_L,wtbr) ;
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
SE = abs(K-repmat(F.',[nM_L 1 nMCMC])).^2 ;
tMSE = (mean(SE,3)) ;
minSE = min(SE,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(SE,[],3) ; ... meandK + std(dK,0,3) ;

CRB = 6/SNR/Ns^3 ;
clf ;
plot(M_L_exp,mindK,':k') ;
plot(M_L_exp,meandK,'.-k','markersize',20) ;
plot(M_L,maxdK,':k') ;
plot(M_L_exp,tMSE,'.r','markersize',35) ;
%plot(M_L,maxSE,':r') ;
%plot(M_L,minSE,':r') ;
plot(M_L_exp,stdK,'ob','markersize',20) ;
set(gca,'xscale','log') ;
set(gca,'yscale','log') ;
axis tight 
grid on



%% 1D UNCERTAINTY ESTIMATION (TWO POLES CLOSELY SPACED IN FREQUENCY)

clc
clear all
%clf ;


Ns = 100 ; % Number of Samples
F0 = 0.0*pi ;
dF = logspace(-4,-1,20)*2*pi ; % Tones (normalized freq.)
FUNC = 'cos' ;
SNR = 1e2 ;
nMCMC = 100 ;
profiler = false ;

F = [F0-dF;F0+dF] ;
U = [100 5] ; % Amplitudes
nF = length(dF) ; % number of tones
switch FUNC
    case 'exp'
        nP = size(F,1) ;
    case 'cos'
        nP = size(F,1)/2 ;
end
t = 0:1:Ns-1 ; % time 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nP ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [] ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = zeros(nF,nP,nMCMC) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for f = 1:nF
    for m = 1:nMCMC
        Amp = ((rand(nP,1)*2-1)+1i*(rand(nP,1)*2-1)) ;
        Amp = Amp./abs(Amp).*U(:) ;
        switch FUNC
            case 'exp'
                Signal = bsxfun(@times,Amp(:),exp(1i*F(:,f)*t)) ;
            case 'cos'
                Signal = bsxfun(@times,Amp(:),cos(F(:,f)*t)) ;
        end
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(sum(Signal+noise,1),EspArgs{:}) ;
        [~,indsort] = sort(real(out.K)) ;
        K(f,:,m) = out.K(indsort) ;
        dK(f,:,m) = abs(out.dK(indsort)).^2 ;
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
mindK = min(abs(dK),[],3) ;
maxdK = max(abs(dK),[],3) ;

SE = abs(bsxfun(@minus,K,F.')).^2 ;
tMSE = mean(SE,3) ;
minSE = min(SE,[],3) ;
maxSE = max(SE,[],3) ;


%plot(dF,mindK,':k') ;
plot(dF,meandK,'.-','markersize',20) ;
%plot(dF,maxdK,':k') ;
%plot(dF,tMSE,'.r','markersize',35) ;
%plot(dF,dF.^-4/Ns^3,'--k','markersize',35) ;
%plot(dF,maxSE,':r') ;
%plot(dF,minSE,':r') ;
%plot(dF,stdK,'ob','markersize',20) ;
plot(dF,dF.^2,'-.k','linewidth',1) ;
set(gca,'xscale','log','yscale','log') ;
grid on




