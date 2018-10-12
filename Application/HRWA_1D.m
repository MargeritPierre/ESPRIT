%% High-Resolution - Wavenumber & frequency analysis

clc
clear all
close all

[file, path] = uigetfile('*.mat','SELECT A FILE FOR HRWA') ;
if path==0; return; end;

DATA = load([path,file]) ;

%% GET INFOS

    % Spatial INFOS
        DATA.nX = 201 ;
        DATA.nY = 3 ;
        DATA.X0 = DATA.X0/1000 ;
        DATA.Y0 = DATA.Y0/1000 ;
        DATA.Lx = DATA.X0(end)-DATA.X0(1) ;
        DATA.Ly = DATA.Y0(end)-DATA.Y0(1) ;
        DATA.dx = DATA.Lx/(DATA.nX-1) ;
        DATA.dy = DATA.Ly/(DATA.nY-1) ;
        DATA.X0 = reshape(DATA.X0,[DATA.nY DATA.nX])-min(DATA.X0(:)) ;
        DATA.Y0 = reshape(DATA.Y0,[DATA.nY DATA.nX])-min(DATA.Y0(:)) ;

    % Time Infos
        if isfield(DATA,'Time')
            DATA.T = DATA.Time ;
            DATA.nT = length(DATA.T) ;
            DATA.nF = DATA.nT ;
            DATA.Lt = DATA.T(end)-DATA.T(1) ;
            DATA.Fe = 1/mean(diff(DATA.T)) ;
            DATA.f = (0:DATA.nF-1)/DATA.nF*DATA.Fe ;
        elseif isfield(DATA,'Freq')
        elseif isfield(DATA,'CorrFreq')
        end
    % Format the displacement
        if isfield(DATA,'Time')
            DATA.U = reshape(permute(DATA.dZ,[2 1]),[DATA.nY DATA.nX DATA.nT]) ;
            DATA.fftU = fft(DATA.U,[],3) ;
        elseif isfield(DATA,'Freq')
        elseif isfield(DATA,'CorrFreq')
        end
        
        
%% MEAN VALUES
    clf ;
    mysubplot(2,1,1) ;
        plot(DATA.T,squeeze(sqrt(mean(mean(DATA.U.^2,1),2))),'k')
        grid on
    mysubplot(2,1,2) ;
        plot(DATA.f,squeeze(sqrt(mean(mean(abs(DATA.fftU).^2,1),2))),'k')
        set(gca,'xscale','log','yscale','log')
        set(gca,'xlim',[1/DATA.nF 1/2]*DATA.Fe)
        grid on
        
        
%% 1D HRWA (wavenumber, harmonic regime)
    % PARAMETERS
        indFmin = 2 ;
        indFmax = find(DATA.f>8e3,1,'first') ; round(DATA.nF/2/1.5) ;
        followSubspace = true ;
        EspArgs = { ...
                    'DIMS_K',[2] ...
                    ,'R0',1:6 ...
                    ,'SOLVER','eig' ...
                    ,'CRITERION','ESTER' ...
                    ,'CRIT_THRS', 1 ...
                    ,'M/L',2/3 ...
                    ,'FUNC','exp' ...
                    ,'SHIFTS',[1]' ...
                    ,'DECIM', [1 1]  ...
                    ,'COMPUTE_U', false  ...
                    ,'COMPUTE_dK', true  ...
                    ,'DEBUG', false  ...
                    }.' ;
        
    % INIT FIGURE
        clf ;
        ax(1) = mysubplot(2,1,1) ;
            ptsK = scatter(NaN,NaN,50,'filled') ;
            %set(gca,'xscale','log')
            set(gca,'yscale','log')
            grid on
            colorbar
        ax(2) = mysubplot(2,1,2) ;
            ptsG = scatter(NaN,NaN,50,'filled') ;
            %set(gca,'xscale','log')
            set(gca,'yscale','log')
            grid on
            colorbar
        colormap(linspecer(100)) ;
        
        
    K = (NaN+1i*NaN)*zeros(DATA.nF,10) ;
    dK = (NaN+1i*NaN)*zeros(DATA.nF,10) ;
    Error = (NaN+1i*NaN)*zeros(DATA.nF,1) ;
    %profile on
    W0 = [] ;
    startTime = tic ;
    t = tic ;
    for ff = indFmin:indFmax % unique(round(linspace(indFmin,indFmax,100))) %
        % Extract Wavevectors
            Signal = DATA.fftU(:,:,ff) ;
            ARGS = {EspArgs{:},'W0',W0} ;
            OUT = ESPRIT(Signal,ARGS{:}) ;
            K(ff,1:size(OUT.K,2)) = OUT.K/DATA.dx ;
            dK(ff,1:size(OUT.dK,2)) = OUT.dK/DATA.dx ;
            Error(ff) = OUT.RelativeError ;
            if followSubspace; W0 = OUT.W ; end
        % Update Figure
            if toc(t)>.5
                t = tic ;
                F = repmat(reshape(DATA.f,[DATA.nF 1]),[1 size(K,2)]) ;
                ptsK.XData = F(:) ;
                ptsK.YData = abs(real(K(:))) ;
                ptsK.CData = log10(sqrt(abs(dK(:)).^2).') ;
                ptsG.XData = F(:) ;
                ptsG.YData = abs(imag(K(:))./real(K(:))) ;
                ptsG.CData = log10(sqrt(abs(dK(:)).^2).') ;
                drawnow ;
            end
    end
    toc(startTime) ;
    %profile viewer
    % Update figure
        valid = abs(imag(K(:))./real(K(:)))<.1 ;
        F = repmat(reshape(DATA.f,[DATA.nF 1]),[1 size(K,2)]) ;
        ptsK.XData = F(valid) ;
        ptsK.YData = abs(real(K(valid))) ;
        ptsK.CData = log10(sqrt(abs(dK(valid)).^2).') ;
        ptsG.XData = F(valid) ;
        ptsG.YData = abs(imag(K(valid))./real(K(valid))) ;
        ptsG.CData = log10(sqrt(abs(dK(valid)).^2).') ;
        axis(ax,'tight')
        drawnow ;

