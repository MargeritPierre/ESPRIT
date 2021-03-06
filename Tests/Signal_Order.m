
%% SIGNAL ORDER ESTIMATION (VARIATION OF THE SNR)

clc
clear all


clf ;
linestyle = '-' ;
Ns = 1000 ; % Number of Samples
nF = 30 ; % number of tones
F = linspace(-.99,.99,nF)*pi*(1+0.01i/Ns*100) ; % Tones (normalized freq.)
NoiseAmp = rand(Ns,1) ; % Noise Amplitudes for non-uniform noise
FUNC = 'exp' ;
SNR = logspace(0,2,20) ;
nMCMC = 10 ;
profiler = false ;

U = [100] ; rand(nF,1) ; % Amplitudes
nSNR = length(SNR) ;
t = 0:1:Ns-1 ; % time 

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 2 ; ...
           'R0' , 1:2*nF ; ...
           'CRITERION' , 'SAMOS' ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , 2/3 ; ...
           'COMPUTE_dK', false ;...
          }' ;

nESTER = zeros(nSNR,nMCMC) ;
nSAMOS = zeros(nSNR,nMCMC) ;
nMDL = zeros(nSNR,nMCMC) ;

wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:nSNR
    for m = 1:nMCMC
        Amp = ((rand(nF,1)*2-1)+1i*(rand(nF,1)*2-1)) ;
        Amp = Amp./abs(Amp).*U(:) ;
        switch FUNC
            case 'exp'
                Signal = bsxfun(@times,Amp(:),exp(1i*F(:)*t)) ;
            case 'cos'
                Signal = bsxfun(@times,Amp(:),cos(F(:)*t)) ;
        end
        Signal = sum(Signal,1) ;
        Rn = rand(1,Ns)*2-1 ;
        In = rand(1,Ns)*2-1 ;
        noise = (Rn+1i*In).*NoiseAmp(:)' ;
        noise = noise*norm(Signal)/norm(noise)/SNR(s) ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        [~,nMDL(s,m)] = max(out.MDL) ;
        [~,nESTER(s,m)] = max(out.ESTER) ;
        [~,nSAMOS(s,m)] = max(out.SAMOS) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/nSNR,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

pMDL = sum(nMDL==nF,2)/nMCMC ;
pESTER = sum(nESTER==nF,2)/nMCMC ;
pSAMOS = sum(nSAMOS==nF,2)/nMCMC ;

plot(SNR,pMDL,linestyle) ;
plot(SNR,pESTER,linestyle) ;
plot(SNR,pSAMOS,linestyle) ;
set(gca,'xscale','log') ;
set(gca,'ylim',[0 1]) ;
lgd = findobj(gcf,'type','legend') ;
if isempty(lgd)
    legend('MDL','ESTER','SAMOS')
else
    strLgd = lgd.String ;
    strLgd{end+1} = 'MDL' ;
    strLgd{end+1} = 'ESTER' ;
    strLgd{end+1} = 'SAMOS' ;
    legend(strLgd) ;
end
grid on


%% SIGNAL ORDER ESTIMATION (TWO POLES CLOSELY SPACED IN FREQUENCY)

clc
clear all


clf ;
Ns = 100 ; % Number of Samples
F0 = 0.0*pi ;
dF = logspace(-4,-1,20)*2*pi ; % Tones (normalized freq.)
FUNC = 'exp' ;
SNR = 1e1 ;
nMCMC = 100 ;
profiler = false ;

F = [F0-dF;F0+dF] ;
U = [100 100] ; % Amplitudes

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
           'R0' , 1:10 ; ...
           'CRITERION' , 'ALL' ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [] ; ...
           'COMPUTE_dK', false ;...
          }' ;

nESTER = zeros(nF,nMCMC) ;
nSAMOS = zeros(nF,nMCMC) ;
nMDL = zeros(nF,nMCMC) ;
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
        [~,nMDL(f,m)] = max(out.MDL) ;
        [~,nESTER(f,m)] = max(out.ESTER) ;
        [~,nSAMOS(f,m)] = max(out.SAMOS) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(f-1)*nMCMC)/nMCMC/nF,wtbr) ;
        end
    end
end
delete(wtbr) ;
if profiler ; profile viewer ; end
drawnow ;

pMDL = sum(nMDL==nP,2)/nMCMC ;
pESTER = sum(nESTER==nP,2)/nMCMC ;
pSAMOS = sum(nSAMOS==nP,2)/nMCMC ;

plot(dF/2/pi,pMDL) ;
plot(dF/2/pi,pESTER) ;
plot(dF/2/pi,pSAMOS) ;
set(gca,'xscale','log') ;
legend('MDL','ESTER','SAMOS')
grid on



