%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SNR)

clc
clear all
%close all
clf ;


Ns = 100 ; % Number of Samples
F = 0.3*pi;%+0.1i ; (linspace(-1,1,1)*.04 + .03)*pi*(1-.10i) ; %rand(1,8)*0.49*pi ; %[.05 -.1 .3 -.02]*pi%] ; % Tones (normalized freq.)];%]%
a = [1+1i] ; % Amplitudes
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
           'DECIM' , [1 2] ; ...
           'FUNC' , FUNC ; 
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [] ; ...
           'COMPUTE_dK', true ;...
           'COMPUTE_U', true ;...
           'COMPUTE_dU', true ;...
          }' ;

K = zeros(nSNR,nF,nMCMC) ;
dK = K ;
Amps = zeros(nSNR,nF,nMCMC) ;
U = zeros(nSNR,nF,nMCMC) ;
dU = zeros(nSNR,nF,nMCMC) ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(1,nF)*2-1)+1i*(rand(1,nF)*2-1)) ;
        Amp = Amp./abs(Amp)*a ;
        Amp = a ; % Comment if you want a random phase
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
        Amps(s,:,m) = Amp ;
        U(s,:,m) = out.U(ind) ;
        dU(s,:,m) = out.dU(ind).^2 ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nSNR,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

% Wavevectors
    % Standard deviation
        stdK = mean(var(K,0,3),2) ;
    % Computed values
        dK = reshape(dK,nSNR,[]) ;
        meandK = mean(abs(dK),2) ;
        mindK = min(abs(dK),[],2) ;
        maxdK = max(abs(dK),[],2) ;
    % Mean Error
        [~,ind] = sort(real(F)) ;
        ME = abs(K-repmat(sort(F(ind)),[nSNR 1 nMCMC])).^2 ;
        ME = reshape(ME,nSNR,[]) ;
        tMSE = (mean(ME,2)) ;
        minSE = min(ME,[],2) ;
        maxSE = max(ME,[],2) ;
    % Figure
        axK = findobj(gcf,'tag','axK') ; if isempty(axK); axK = mysubplot(1,2,1) ; end
            title 'Uncertainty on K'
            plot(SNR,mindK,':k') ;
            plot(SNR,meandK,'.-k','markersize',20) ;
            plot(SNR,maxdK,':k') ;
            plot(SNR,tMSE,'.r','markersize',35) ;
            plot(SNR,maxSE,':r') ;
            plot(SNR,minSE,':r') ;
            plot(SNR,stdK,'ob','markersize',20) ;
            set(gca,'xscale','log','yscale','log') ;
            grid on

% Amplitudes
    % Standard deviation
        stdU = mean(var(U,0,3),2) ;
    % Computed values
        dU = reshape(dU,nSNR,[]) ;
        meandU = mean(abs(dU),2) ;
        mindU = min(abs(dU),[],2) ;
        maxdU = max(abs(dU),[],2) ;
    % Mean Error
        MEU = abs(U-Amps).^2 ;
        MEU = reshape(MEU,nSNR,[]) ;
        tMSEU = (mean(MEU,2)) ;
        minSEU = min(MEU,[],2) ;
        maxSEU = max(MEU,[],2) ;
    % Figure
        axU = findobj(gcf,'tag','axU') ; if isempty(axU); axU = mysubplot(1,2,2) ; end
            title 'Uncertainty on U'
            plot(SNR,mindU,':k') ;
            plot(SNR,meandU,'.-k','markersize',20) ;
            plot(SNR,maxdU,':k') ;
            plot(SNR,tMSEU,'.r','markersize',35) ;
            plot(SNR,maxSEU,':r') ;
            plot(SNR,minSEU,':r') ;
            plot(SNR,stdU,'ob','markersize',20) ;
            set(gca,'xscale','log','yscale','log') ;
            grid on



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SIGNAL SIZE)

clc
clear all
%close all
clf ;


Ns = round(logspace(log10(20),log10(100),20)) ; % Number of Samples
F = [.3]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
FUNC = 'cos' ;
SNR = 1e2 ;
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
           'M/L' , [] ; ...
           'COMPUTE_dK', true ;...
           'COMPUTE_U', true ;...
           'COMPUTE_dU', true ;...
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
ME = abs(K-repmat(F.',[nNs 1 nMCMC])).^2 ;
tMSE = (mean(ME,3)) ;
minSE = min(ME,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(ME,[],3) ; ... meandK + std(dK,0,3) ;


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

clf 

Ns = 100 ; % Number of Samples
F = logspace(log10(0.1),log10(Ns/10),10)*pi/Ns*(1+.00i) ; % Tones (normalized freq.)
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
           'M/L' , [] ; ...
           'COMPUTE_dK', true ;...
           'COMPUTE_U', true ;...
           'COMPUTE_dU', true ;...
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
        Amp = U ;
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
ME = abs(K-repmat(F.',[1 1 nMCMC])).^2 ;
ME = reshape(ME,nF,[]) ;
tMSE = (mean(ME,2)) ;
minSE = min(ME,[],2) ; ... meandK - std(dK,0,3) ;
maxSE = max(ME,[],2) ; ... meandK + std(dK,0,3) ;

plot(F/2/pi*Ns,mindK,':k') ;
plot(F/2/pi*Ns,meandK,'.-k','markersize',20) ;
plot(F/2/pi*Ns,maxdK,':k') ;
plot(F/2/pi*Ns,tMSE,'.r','markersize',35) ;
plot(F/2/pi*Ns,maxSE,':r') ;
plot(F/2/pi*Ns,minSE,':r') ;
plot(F/2/pi*Ns,stdK,'ob','markersize',20) ;
%plot(F/pi,(F/pi).^2,'-.r','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;
grid on



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SNR, Multiple Snapshot)

clc
clear all
%close all


Ns = 100 ; % Number of Sensors
Nsnap = 10 ; % Number of Snapshots
FUNC = 'exp' ;
F = [-.023-.01i]*pi ; % Tones (normalized freq.)
U = [1] ; % Amplitudes
SNR = logspace(-1,4,10) ;
nMCMC = 100 ;
profiler = false ;

t = 0:1:Ns-1 ; % time 
nF = length(F) ; % number of tones
nSNR = length(SNR) ; % number of SNR levels
if strcmp(FUNC,'cos'); F = F.*sign(real(F)) ; end

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , nF ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'TLS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , 1/2 ; ... 
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
        switch FUNC
            case 'exp'
                Signal = Amp*sum(exp(1i*F.'*t),1) ; 
            case 'cos'
                Signal = Amp*sum(cos(1i*F.'*t),1) ; 
        end
        Signal = Amp*sum(exp(1i*F.'*t),1) ; 
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        K(s,:,m) = out.K ;
        dK(s,:,m) = out.dK.^2 ;
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
MSE = abs(K-repmat(F.',[nSNR 1 nMCMC])).^2 ;
tMSE = mean(MSE,3) ;
minSE = min(MSE,[],3) ;
maxSE = max(MSE,[],3) ;


clf ;
plot(SNR,mindK,':k') ;
plot(SNR,meandK,'.-k','markersize',20) ;
plot(SNR,maxdK,':k') ;
plot(SNR,tMSE,'.r','markersize',35) ;
plot(SNR,maxSE,':r') ;
plot(SNR,minSE,':r') ;
plot(SNR,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;



%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE NUMBER OF SNAPSHOTS)

clc
clear all
%close all


Ns = 50 ; % Number of Sensors
Nsnap = unique(round(logspace(log10(5),log10(50),10))) ; % Number of Snapshots
F = [-.000023-.1i*0]*pi ; % Tones (normalized freq.)
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
           'FIT' , 'TLS' ; ...
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
        dK(s,:,m) = out.dK.^2 ;
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

ME = abs(K-repmat(F.',[nNsnap 1 nMCMC])).^2 ;
tMSE = (mean(ME,3)) ;
minSE = min(ME,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(ME,[],3) ; ... meandK + std(dK,0,3) ;

clf ;
plot(Nsnap,mindK,':k') ;
plot(Nsnap,meandK,'.-k','markersize',20) ;
plot(Nsnap,maxdK,':k') ;
plot(Nsnap,tMSE,'.r','markersize',35) ;
plot(Nsnap,maxSE,':r') ;
plot(Nsnap,minSE,':r') ;
plot(Nsnap,stdK,'ob','markersize',20) ;
set(gca,'xscale','log','yscale','log') ;








%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE SPATIAL SMOOTHING)

clc
clear all
%close all


clf ;
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
        dK(k,:,m) = out.dK^2 ;
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
ME = abs(K-repmat(F.',[nM_L 1 nMCMC])).^2 ;
tMSE = (mean(ME,3)) ;
minSE = min(ME,[],3) ; ... meandK - std(dK,0,3) ;
maxSE = max(ME,[],3) ; ... meandK + std(dK,0,3) ;

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

ME = abs(bsxfun(@minus,K,F.')).^2 ;
tMSE = mean(ME,3) ;
minSE = min(ME,[],3) ;
maxSE = max(ME,[],3) ;


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




%% 1D UNCERTAINTY ESTIMATION (VARIATION OF THE NUMBER OF COMPONENTS)

clc
clear all


%clf ;
linestyle = '-' ;
decim = 1 ;
Ns = 100*decim ; % Number of Samples
nT = 1:15 ; unique(ceil(linspace(1,1/4*Ns,10))) ; % number of tones
Fmax = .8/decim ; 
Fmin = 0.05/decim ; 0.05 ;  
eta = 0.01 ;
SNR = 1e2 ;
nMCMC = 100 ;
NoiseAmp = 1 ; rand(Ns,1) ; % Noise Amplitudes for non-uniform noise
profiler = false ;

U = [100] ; % Amplitudes
t = 0:1:Ns-1 ; % time 
nTmax = length(nT) ;

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'CRITERION' , 'MDL' ; ...
           'DECIM' , [1 decim] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , 1/2 ; ...
           'COMPUTE_dK', true ;...
          }' ;

K = (NaN+1i*NaN)*zeros(nTmax,nTmax,nMCMC) ;
dK = (NaN)*zeros(nTmax,nTmax,nMCMC) ; 
ME = (NaN)*zeros(nTmax,nTmax,nMCMC) ;  
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for ff = 1:nTmax
    nF = nT(ff) ;
    F = linspace(Fmin,Fmax,nF)*pi*(1+1i*eta/Ns*100) ; % Tones (normalized freq.)
    for m = 1:nMCMC
        Amp = ((rand(nF,1)*2-1)+1i*(rand(nF,1)*2-1)) ;
        Amp = Amp./abs(Amp).*U(:) ;
        Signal = bsxfun(@times,Amp(:),exp(1i*F(:)*t)) ;
        Signal = sum(Signal,1) ;
        Rn = rand(1,Ns)*2-1 ;
        In = rand(1,Ns)*2-1 ;
        noise = (Rn+1i*In).*NoiseAmp(:)' ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        ESPARGS = {EspArgs{:},'R0',nF} ;
        out = ESPRIT(Signal+noise,ESPARGS{:}) ;
        [K(ff,1:nF,m),indsort] = sort(out.K) ;
        dK(ff,1:nF,m) = out.dK(indsort).^2 ;
        ME(ff,1:nF,m) = abs(K(ff,1:nF,m)-F).^2 ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(ff-1)*nMCMC)/nMCMC/nTmax,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;


meandK = mean(dK,3) ;
MSE = mean(ME,3) ;

notNaN = ~isnan(K(:,:,1)) ;
totalMeandK = meandK ; totalMeandK(isnan(totalMeandK)) = 0 ; totalMeandK = sum(bsxfun(@times,totalMeandK,notNaN),2)./sum(notNaN,2) ;
totalMSE = MSE ; totalMSE(isnan(totalMSE)) = 0 ; totalMSE = sum(bsxfun(@times,totalMSE,notNaN),2)./sum(notNaN,2) ;

%plot(1:nFmax,meandK)
%plot(1:nFmax,MSE,'ob')
plot(1:nTmax,totalMeandK,'-.')
set(gca,'colororderindex',get(gca,'colororderindex')-1) ;
plot(1:nTmax,totalMSE,'.','markersize',35)
set(gca,'yscale','log')
grid on


