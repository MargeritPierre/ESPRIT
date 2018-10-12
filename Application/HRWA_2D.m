clc
clear all
close all

[file, path] = uigetfile('*.mat','SELECT A FILE FOR HRWA') ;
if path==0; return; end;

DATA = load([path,file]) ;

DATA.U3 = double(DATA.U3) ;
[DATA.nY,DATA.nX,DATA.nF] = size(DATA.U3) ;
DATA.dx = DATA.XX(2,2)-DATA.XX(1,1) ;
DATA.dy = DATA.YY(2,2)-DATA.YY(1,1) ;
DATA.kx = indfftshift((0:DATA.nX-1)/DATA.nX*2*pi/DATA.dx) ;
DATA.ky = indfftshift((0:DATA.nY-1)/DATA.nY*2*pi/DATA.dy) ;
[DATA.KX,DATA.KY] = meshgrid(DATA.kx,DATA.ky) ;
DATA.fftU3 = fftshift(fftshift(fft2(DATA.U3),1),2) ;

DATA


%% FIRST EXTRACTION
    % PARAMETERS
        indFmin = 1 ;
        indFmax = round(DATA.nF/2/1.2) ;
        DECIM = [1 1]*1 ;
        followSubspace = false ;
        EspArgs = { ...
                    'DIMS_K',[1 2] ...
                    ,'R0',1:10 ...
                    ,'SOLVER','eig' ...
                    ,'CRITERION','ESTER' ...
                    ,'CRIT_THRS', 1 ...
                    ,'M/L',2/3 ...
                    ,'FUNC',repmat('exp',[2 1]) ...
                    ,'SHIFTS',eye(2) ...
                    ,'DECIM', DECIM  ...
                    ,'COMPUTE_U', false  ...
                    ,'COMPUTE_dK', true  ...
                    ,'DEBUG', true  ...
                    }.' ;
        
    % INIT FIGURE
        clf ;
        srf = surf(DATA.KX,DATA.KY,DATA.KY*0 ...
                ...,real(DATA.U3(:,:,1072))...
                ,log10(abs(DATA.fftU3(:,:,784)))...
                ,'facecolor','interp'...
                ,'edgecolor','none'...
                ,'facealpha',0.5...
                ) ;
        ptsK = scatter3(NaN,NaN,NaN,50,'filled') ; %plot3(NaN,NaN,NaN,'.k','markersize',20) ;
        set(gca,'view',[32 32]) ;
        %axis tight
        set(gca,'xlim',DATA.kx([1,end])/DECIM(2))
        set(gca,'ylim',DATA.ky([1,end])/DECIM(1))
        set(gca,'zlim',[DATA.f(indFmin) DATA.f(indFmax)])
        colorbar
        colormap(linspecer(100)) ;
        % myaxisequal('xy')
        
        
    K = (NaN+1i*NaN)*zeros(2,100,DATA.nF) ;
    dK = (NaN+1i*NaN)*zeros(2,100,DATA.nF) ;
    %profile on
    W0 = [] ;
    startTime = tic ;
    for ff = indFmin:indFmax % unique(round(linspace(indFmin,indFmax,100))) %
        % Extract Wavevectors
            Signal = DATA.U3(:,:,ff) ;
            ARGS = {EspArgs{:},'W0',W0} ;
            OUT = ESPRIT(Signal,ARGS{:}) ;
            K(:,1:size(OUT.K,2),ff) = diag(1./[DATA.dy,DATA.dx])*OUT.K ;
            dK(:,1:size(OUT.dK,2),ff) = diag(1./[DATA.dy,DATA.dx])*OUT.dK ;
            if followSubspace; W0 = OUT.W ; end
        % Update Figure
            F = repmat(reshape(DATA.f,[1 DATA.nF]),[size(K(1,:,1),2) 1]) ;
            srf.CData = log10(abs(DATA.fftU3(:,:,ff))) ;
            srf.ZData = srf.ZData*0 + DATA.f(ff) ;
            ptsK.XData = real(K(2,:)).' ;
            ptsK.YData = real(K(1,:)).' ;
            ptsK.ZData = F(:) ;
            ptsK.CData = log10(sqrt(sum(abs(dK(:,:)).^2,1)).') ;
            drawnow ;
    end
    toc(startTime) ;
    %profile viewer
    pause(.5)
    delete(srf) ;
    myaxisequal('xy',1)
    axis tight
    %set(gca,'xscale','log','yscale','log','zscale','log')
    
    
    
%% EQUIVALENT BENDING STIFFNESS

    Dmax = 210e9*(10e-3)^2/7800 ;

    kx = real(K(2,:)) ;
    ky = real(K(1,:)) ;
    dkx = abs(dK(2,:))/2 ;
    dky = abs(dK(1,:))/2 ;
    f = F(:).' ;
    
    k = sqrt(kx.^2+ky.^2) ;
    dk = (abs(dkx.*kx)+abs(dky.*ky))./k ;
    theta = angle(kx+1i*ky) ;
    dtheta = 2*dk./k ;
    
    D = (2*pi*f).^2./k.^4 ;
    
    % VALID POINTS
        valid = ~isnan(D) ;
        valid = valid & abs(D)<Dmax ;
    
    
    clf
    scat = scatter3(D(valid).*cos(theta(valid)),D(valid).*sin(theta(valid)),f(valid),25,'filled') ;
    scat.CData = log10(dk(valid)) ;
    colorbar
    axis tight ;
    polargrid ;
    myaxisequal('xy',2) ;
    axis off
    
    
%% KIRCHHOFF THEORY
    
    IDENTIF = ... 'polar' ...
               'xyz' ...
                ;
    FIT = ... 'LS' ...
           'TLS' ...
          ... 'GLS' ...
            ;
    
    % Linear system A.x = b
        % Matrix A
            switch IDENTIF
                case 'xyz'
                    c = cos(theta(:)) ;
                    s = sin(theta(:)) ;
                    A = [   c.^4 ...
                            s.^4 ...
                            2*c.^2.*s.^2 ...
                            4*c.^3.*s ...
                            4*c.*s.^3
                        ] ;
                    c = abs(c) ;
                    s = abs(s) ;
                    dA = [  4*s.*c.^3.*dtheta(:) ...
                            4*c.*s.^3.*dtheta(:) ...
                            4*(c.*s.^3 + c.^3.*s).*dtheta(:) ...
                            4*(3*c.^2.*s.^2 + c.^4).*dtheta(:) ...
                            4*(3*c.^2.*s.^2 + s.^4).*dtheta(:)
                         ] ;
                case 'polar'
                    A = [   ones(size(k(:))) ...
                            exp(2i*theta(:)) ...
                            exp(-2i*theta(:)) ...
                            exp(4i*theta(:)) ...
                            exp(-4i*theta(:))
                        ] ;
                    dA = [  zeros(size(k(:))) ...
                            2i.*dtheta(:) ...
                            2i.*dtheta(:) ...
                            4i.*dtheta(:) ...
                            4i.*dtheta(:)
                         ] ;
            end
        % Vector b
            b = (2*pi*f(:)).^2./k(:).^4 ;
            db = 4*(2*pi*f(:)).*dk(:).^3./k(:).^8 ;
        % Cull invalid points
            A = A(valid,:) ;
            b = b(valid) ;
            dA = dA(valid,:) ;
            db = db(valid) ;
            
    % Linear regression
        switch FIT
            case 'LS'
                Di = A\b ;
            case 'TLS'
                [Di,w] = eig([A b]'*[A b],'vector') ;
                [~,indsort] = sort(w,'descend') ;
                Di = Di(:,indsort) ;
                Di = -Di(1:5,end)/Di(end,end) ;
            case 'GLS'
                w = 1./dk(valid) ;
                A = diag(w)*A ;
                b = diag(w)*b ;
                Di = A\b ;
        end
        Di
        dDi = abs(A\(db + dA*abs(Di)))
        dDi = abs(A\(A*Di-b))
        dDi = sqrt(var(b)*diag(pinv(A)*pinv(A)'))
        
    % Conversion and model reconstruction
        phi = linspace(0,2*pi,200) ;
        switch IDENTIF
            case 'polar'
                T = real(Di(1)) ;
                R0 = 2*real(sqrt(Di(4)*Di(5))) ;
                R1 = 2*real(sqrt(Di(2)*Di(3))) ;
                phi0 = real(-1i/4*log(Di(5)/Di(4))) ;
                phi1 = real(-1i/2*log(Di(3)/Di(2))) ;
                dT = real(dDi(1)) ;
                dR0 = real(dDi(4)*Di(5)+Di(4)*dDi(5))/R0 ;
                dR1 = 2*real(sqrt(Di(2)*Di(3))) ;
                dphi0 = real(-1i/4*log(Di(5)/Di(4))) ;
                dphi1 = real(-1i/2*log(Di(3)/Di(2))) ;
                Dm = T+R0*cos(4*(phi-phi0))+R1*cos(2*(phi-phi1));
            case 'xyz'
                Dm = Di(1)*cos(phi).^4 ...
                        + Di(2)*sin(phi).^4 ...
                        + 2*Di(3)*(cos(phi).^2.*sin(phi).^2) ...
                        + 4*Di(4)*(cos(phi).^3.*sin(phi).^1) ...
                        + 4*Di(5)*(cos(phi).^1.*sin(phi).^3) ...
                    ;
                dDm = dDi(1)*cos(phi).^4 ...
                        + dDi(2)*sin(phi).^4 ...
                        + 2*dDi(3)*(cos(phi).^2.*sin(phi).^2) ...
                        + abs(4*dDi(4)*(cos(phi).^3.*sin(phi).^1)) ...
                        + abs(4*dDi(5)*(cos(phi).^1.*sin(phi).^3)) ...
                    ; 
        end
        
    % PLOT
        delete(findobj(gca,'tag','model')) ;
        plot3(...
                Dm(:).*cos(phi(:)) ...
                ,Dm(:).*sin(phi(:)) ...
                ,repmat([min(f(valid)) max(f(valid))],[length(Dm) 1]) ...
                ,'-k' ...
                ,'linewidth',2 ...
                ,'tag','model' ...
            ) ;
        plot3(...
                (Dm(:)+dDm(:)).*cos(phi(:)) ...
                ,(Dm(:)+dDm(:)).*sin(phi(:)) ...
                ,repmat([min(f(valid)) max(f(valid))],[length(Dm) 1]) ...
                ,'-.k' ...
                ,'linewidth',1 ...
                ,'tag','model' ...
            ) ;
        plot3(...
                (Dm(:)-dDm(:)).*cos(phi(:)) ...
                ,(Dm(:)-dDm(:)).*sin(phi(:)) ...
                ,repmat([min(f(valid)) max(f(valid))],[length(Dm) 1]) ...
                ,'-.k' ...
                ,'linewidth',1 ...
                ,'tag','model' ...
            ) ;
        













    
    
    
    
    
