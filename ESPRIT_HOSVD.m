%% HIGHER-ORDER SVD-BASED ESPRIT
%% OPEN A FILE

clc
clear all
close all

[file,path] = uigetfile('.mat') ;
if file==0; return; end

DATA = load([path,file]) ;

DATA.dx = DATA.XX(1,2) - DATA.XX(1,1) ;
DATA.dy = DATA.YY(2,1) - DATA.YY(1,1) ;
[DATA.nY,DATA.nX,DATA.nF] = size(DATA.U3) ;

DATA

%% Plot mean values
    close all
    fig = figure('windowstyle','docked');
    DATA.meanU3 = mean(abs(reshape(DATA.U3,[DATA.nY*DATA.nX DATA.nF])),1) ;
    axF = axes() ;
    plot(DATA.meanU3(1:floor(end/2)),'k')
    set(gca,'xscale','log','yscale','log')
    axis tight
    grid on
    crs = plot([1 1],get(gca,'ylim'),'r') ;
    ttl = title('1') ;
    ax = axes('position',[0.1 0.1 .4 .4]) ;
    uistack(ax,'bottom') ;
    axis off
    srf = surf(DATA.XX,DATA.YY,real(DATA.U3(:,:,1))) ;
    shading interp ;
    axis equal
    axis tight
    mousepos = @()max(min(axF.XLim),min(round(axF.CurrentPoint(1)),max(axF.XLim))) ;
    fig.WindowButtonMotionFcn = @(src,evt){set(...
                                                crs,'XData',mousepos()*[1 1]...
                                                );...
                                            set(...
                                                ttl,'string',['i = ',num2str(mousepos()),', f = ',num2str(DATA.f(mousepos()),4),' Hz']...
                                                );...
                                            set(...
                                                srf,'ZData',real(DATA.U3(:,:,mousepos()))...
                                                );...
                                            };


%% SVD DECOMPOSITION VISU
    close all
    ff = 384 ;
    [U,S,V] = svd(DATA.U3(:,:,ff)) ;
    s = diag(S) ;
    %criterion = cumsum(s)/sum(s) ;
    criterion = [s(1:end-1)-s(2:end);0]./s ;
    %criterion = diag(s)/s(1) ;

    fig = figure('windowstyle','docked');
    axF = axes() ;
    plot(criterion,'k')
    axis tight
    %set(gca,'xscale','log')
    %set(gca,'yscale','log')
    %set(gca,'ylim',[0 1]) ;
    grid on
    crs = plot([1 1],get(gca,'ylim'),'r') ;
    ttl = title('1') ;
    ax = axes('position',[0.025 0.025 .45 .45]) ;
    %ax.View = [30 30]
    uistack(ax,'bottom') ;
    axis off
    srf = surf(DATA.XX,DATA.YY,real(DATA.U3(:,:,ff))) ;
    shading interp ;
    axis equal
    axis tight
    axErr = axes('position',[0.525 0.525 .45 .45]) ;
    uistack(axErr,'bottom') ;
    axis off
    srfErr = surf(DATA.XX,DATA.YY,real(DATA.U3(:,:,ff))) ;
    shading interp ;
    axis equal
    axis tight
    mousepos = @()max(min(axF.XLim),min(round(axF.CurrentPoint(1)),max(axF.XLim))) ;
    constru = @(r)U(:,1:r)*S(1:r,1:r)*V(:,1:r)' ;
    fig.WindowButtonMotionFcn = @(src,evt){set(...
                                                crs,'XData',mousepos()*[1 1]...
                                                );...
                                            set(...
                                                ttl,'string',['i = ',num2str(mousepos()),', f = ',num2str(DATA.f(mousepos()),4),' Hz']...
                                                );...
                                            set(...
                                                srf,'ZData',real(constru(mousepos()))...
                                                );...
                                            set(...
                                                srfErr,'ZData',abs(DATA.U3(:,:,ff)-constru(mousepos()))...
                                                );...
                                            };
                                        
                                        
%% SVD AND DERIVATION

    indF = 10:440 ; 20:300 ; 50:1:1500 ;
    Rmin = 1 ; 
    Rmax = round(min(DATA.nX,DATA.nY)*1/2) ;3/4
    mdl_thrs = 1 ; .65 ;
    
    Mx = ceil(DATA.nX/2) ;
    My = ceil(DATA.nY/2) ;
    Kx = DATA.nX-Mx+1 ;
    Ky = DATA.nY-My+1 ;
    
    nIndF = length(indF) ;
    wtbr = waitbar(0) ;
    imdl = zeros(1,nIndF)*NaN ;
    mdl = zeros(1,nIndF)*NaN ;
    t = tic ;
    KX = [] ;
    KY = [] ;
    F = [] ;
    for ff = 1:1:nIndF
        % Signal Retrieving
            sig = DATA.U3(:,:,indF(ff)) ;
        % Spatial Smoothing
            covar = sig(1:Ky,1:Kx) ;
            for ii = 2:My
                for jj = 2:Mx
                    covar = covar+sig((1:Ky)+ii-1,(1:Kx)+jj-1) ;
                end
            end
        % SVD Decomposition
            [U,S,V] = svd(covar) ;
            s = diag(S) ;
        % MDL CRITERION
            MDL = [[s(1:end-1)-s(2:end);0]./s] ;
            MDL = MDL(Rmin:Rmax) ;
            maxMDL = max(MDL) ;
            imdl(ff) = find(MDL>=mdl_thrs*maxMDL,1,'last')+Rmin-1 ;
            mdl(ff) = MDL(imdl(ff)-Rmin+1) ;
        % TRUNCATION
            R = imdl(ff) ;
            U = U(:,1:R) ;
            V = V(:,1:R) ;
            S = S(1:R,1:R) ;
        % SHIFTED SIGNALS
            WupX = kron(U,ones(Kx-1,1)).*kron(ones(Ky,1),V(1:end-1,:)) ;
            WdwnX = kron(U,ones(Kx-1,1)).*kron(ones(Ky,1),V(2:end,:)) ;
            WupY = kron(U(1:end-1,:),ones(Kx,1)).*kron(ones(Ky-1,1),V) ;
            WdwnY = kron(U(2:end,:),ones(Kx,1)).*kron(ones(Ky-1,1),V) ;
        % SPECTRAL MATRICES
            Fx = WupX\WdwnX ; % V(1:end-1,:)\V(2:end,:) ; % WupX\WdwnX ; %
            Fy = WupY\WdwnY ; % U(1:end-1,:)\U(2:end,:) ; % WupY\WdwnY ; %
        % POLAR MATRICES
           [T,~] = eig(Fx+1.1*Fy) ;% T = eye(R) ; %
            Zx = T\Fx*T ;
            Zy = T\Fy*T ;
        % WAVENUMBERS
            KX(end+(1:R)) = -1i*log(diag(Zx))/DATA.dx ;
            KY(end+(1:R)) = -1i*log(diag(Zy))/DATA.dy ;
            F(end+(1:R)) = DATA.f(ff) ;
        % WAITBAR
            wtbr = waitbar(ff/nIndF,wtbr) ;
    end
    toc(t) 
    delete(wtbr) ;
    
    clf('reset') ;
    %plot(indF,imdl)
    plot3(real(KX),real(KY),F,'.k') ;
    myaxisequal('xy')
    
    
%% SPATIAL DESCRIPTION

D = zeros(nYc,nXc,5) ;
for xx = 1:nXc
    for yy = 1:nYc
        A = [dw_dx4(yy,xx,:)...
            dw_dy4(yy,xx,:)...
            2*dw_dx2dy2(yy,xx,:)...
            4*dw_dx3dy(yy,xx,:)...
            4*dw_dxdy3(yy,xx,:)...
            ] ;
        b = (2*pi*DATA.f(indF).').^2.*squeeze(w(yy,xx,:)) ;
        A = squeeze(A).' ;
        weights = 1;%./b.^5;%./b ;
        D(yy,xx,:) = (diag(weights)*A)\(diag(weights)*b) ;
    end
end


clf('reset') ;
ax = mysubplot(5,1,1) ;
    srf = surf(DATA.XX(1:nYc,1:nXc),DATA.YY(1:nYc,1:nXc),real(D(:,:,1))) ;
    ttl = title('$B_{11}$') ;
    colorbar ;
ax(end+1) = mysubplot(5,1,2) ;
    srf(end+1) = surf(DATA.XX(1:nYc,1:nXc),DATA.YY(1:nYc,1:nXc),real(D(:,:,2))) ;
    ttl(end+1) = title('$B_{22}$') ;
    colorbar ;
ax(end+1) = mysubplot(5,1,3) ;
    srf(end+1) = surf(DATA.XX(1:nYc,1:nXc),DATA.YY(1:nYc,1:nXc),real(D(:,:,3))) ;
    ttl(end+1) = title('$B_{12}+ 2 B_{66}$') ;
    colorbar ;
ax(end+1) = mysubplot(5,1,4) ;
    srf(end+1) = surf(DATA.XX(1:nYc,1:nXc),DATA.YY(1:nYc,1:nXc),real(D(:,:,4))) ;
    ttl(end+1) = title('$B_{16}$') ;
    colorbar ;
ax(end+1) = mysubplot(5,1,5) ;
    srf(end+1) = surf(DATA.XX(1:nYc,1:nXc),DATA.YY(1:nYc,1:nXc),real(D(:,:,5))) ;
    ttl(end+1) = title('$B_{26}$') ;
    colorbar ;
set(srf,'facecolor','interp','edgecolor','none') ;
set(ax,'xtick',[],'ytick',[],'ztick',[])
set(ax,'box','on')
axis(ax,'equal')
axis(ax,'tight')


%% DELETE EXTREMUMS

for i=1:5
   data = real(D(:,:,i)) ;
   med = median(data(:)) ;
   ec = var(sqrt(data(:)),0,1) ;
   valid = abs(data-med)<2*ec ;
   caxis(ax(i),med+5*ec*[-1 1]) ;
end


%% FREQUENCY DESCRIPTION

D = zeros(nIndF,5) ;
for ff = 1:nIndF
        A = cat(3,...
            dw_dx4(:,:,ff)...
            ,dw_dy4(:,:,ff)...
            ,2*dw_dx2dy2(:,:,ff)...
            ,4*dw_dx3dy(:,:,ff)...
            ,4*dw_dxdy3(:,:,ff)...
            ) ;
        b = (2*pi*DATA.f(ff).').^2.*w(:,:,ff) ;
        A = reshape(A,[],5) ; b = b(:) ;
        weights = 1 ; %b ;
        D(ff,:) = (diag(weights)*A)\(diag(weights)*b) ;
end

clf('reset') ;
mysubplot(2,1,1) ;
    plot(DATA.f(indF),real(D))
mysubplot(2,1,2) ;
    plot(DATA.f(indF),imag(D(:,1:2))./real(D(:,1:2)))












