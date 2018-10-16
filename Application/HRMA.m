%% PROCESSING THE MEASUREMENTS OF THE CANTILEVER S BEAM

clc
clear all

[file,path] = uigetfile('*.mat','SELECT A TRANSIENT MEASUREMENT FILE FOR HRMA') ;
if path==0; return; end

DATA = load([path,file]) ;

%% GET INFOS

    % Spatial INFOS
        DATA.nX = 201 ;
        DATA.nY = 3 ;
        DATA.X0 = DATA.X0/1000 ;
        DATA.Y0 = DATA.Y0/1000 ;
        DATA.Z0 = DATA.Z0/1000 ;
        DATA.Lx = DATA.X0(end)-DATA.X0(1) ;
        DATA.Ly = DATA.Y0(end)-DATA.Y0(1) ;
        DATA.Lz = DATA.Z0(end)-DATA.Z0(1) ;
        DATA.dx = DATA.Lx/(DATA.nX-1) ;
        DATA.dy = DATA.Ly/(DATA.nY-1) ;
%         DATA.X0 = reshape(DATA.X0,[DATA.nY DATA.nX])-min(DATA.X0(:)) ;
%         DATA.Y0 = reshape(DATA.Y0,[DATA.nY DATA.nX])-min(DATA.Y0(:)) ;
%         DATA.Z0 = reshape(DATA.Z0,[DATA.nY DATA.nX])-min(DATA.Z0(:)) ;

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
        DATA.meanNRJ = squeeze(sqrt(mean(mean(DATA.U.^2,1),2))) ;
        plot(DATA.T,DATA.meanNRJ,'k')
        DATA.upEdg = find(DATA.meanNRJ>max(DATA.meanNRJ)*0.1,1,'first') ;
        plot(DATA.T(DATA.upEdg)*[1 1],get(gca,'ylim'),':r')
        grid on
    mysubplot(2,1,2) ;
        plot(DATA.f,squeeze(sqrt(mean(mean(abs(DATA.fftU).^2,1),2))),'k')
        set(gca,'xscale','log','yscale','log')
        set(gca,'xlim',[1/DATA.nF 1/2]*DATA.Fe)
        grid on

%% TIME PREVIEW

    % Params
        indT = DATA.upEdg-10:DATA.nT ;
        scaleFactor = .15 ; .06 ;


    NORM = norm([DATA.Lx DATA.Ly DATA.Lz]) ;
    maxAmp = NORM*scaleFactor ;
    dPx = DATA.dX ;
    dPy = DATA.dY ;
    dPz = DATA.dZ ;
    U = sqrt(dPx.^2+dPy.^2+dPz.^2) ;
    maxU = max(U(:)) ;
    dXn = dPx./maxU.*maxAmp ;
    dYn = dPy./maxU.*maxAmp ;
    dZn = dPz./maxU.*maxAmp ;
    
    ax = gca ; 
        ax.Position = [0 0 1 1] ;
        delete(ax.Children) ;
        pl0 = plot3(DATA.X0,DATA.Y0,DATA.Z0,'-.k','linewidth',.5);
        pl = plot3(DATA.X0,DATA.Y0,DATA.Z0,'.r','markersize',10);
        axis equal ;
        ax.XLim = [min(DATA.X0(:)) max(DATA.X0(:))] + [-1 1]*maxAmp*1.1 ;
        ax.YLim = [min(DATA.Y0(:)) max(DATA.Y0(:))] + [-1 1]*maxAmp*1.1 ;
        ax.ZLim = [min(DATA.Z0(:)) max(DATA.Z0(:))] + [-1 1]*maxAmp*1.1 ;
        %ax.View = [145 40] ;
    
    for tt = indT
        pl.XData = DATA.X0 + dXn(tt,:).' ;
        pl.YData = DATA.Y0 + dYn(tt,:).' ;
        pl.ZData = DATA.Z0 + dZn(tt,:).' ;
        drawnow ;
        pause(0.02) ;
    end
        
%% MODAL ANALYSIS

    % Create Signal
        indT =  DATA.upEdg+100+(1:1000) ;
        Signal = [DATA.dZ].' ;
        Signal = Signal(:,indT) ;
        %Signal = Signal-mean(Signal(:)) ;

    % Esprit Parameters

    % APPLY ESPRIT
        OUT = ESPRIT(Signal ...
                        ,'DIMS_K', 2 ...
                        ,'SOLVER', 'eig' ...
                        ,'CRITERION', 'ESTER' ...
                        ,'CRIT_THRS', 1 ...
                        ,'COMPUTE_U', true ...
                        ,'SIGNAL_MODEL', true ...
                        ,'COMPUTE_dK', false ...
                        ,'FUNC', 'exp' ...
                        ,'R0', 2:2:100 ...
                        ,'M/L', 2/3 ...
                        ,'W0', [] ...
                        ,'FIT', 'TLS' ...
                        ,'DECIM', [1 10] ...
                        ,'SHIFTS', [1]' ...
                        ,'DEBUG', true ...
                        ,'STABILDIAG', true ...
                        ,'MAC',.1 ...
                        ) 
        OUT.K = OUT.K*DATA.Fe/2/pi ;
        OUT.dK = OUT.dK*DATA.Fe/2/pi ;
                    
    % RESULTS
        clf ;
            % Measurements
                DATA.meanfftU = squeeze(sqrt(mean(mean(abs(DATA.fftU).^2,1),2))) ;
                plot(DATA.f,DATA.meanfftU,'k','linewidth',1)
            % Signal Model
                fModel = (0:length(indT)-1)/length(indT)*DATA.Fe ;
                plot(fModel,sqrt(mean(abs(fft(OUT.SignalModel,[],2).^2),1)),':r')
            % Formatting
                set(gca,'xscale','log','yscale','log')
                set(gca,'xlim',[1/DATA.nF 1/2]*DATA.Fe)
                grid on
            % Identified Parameters
                set(gca,'colororderindex',1) ;
                set(gca,'colororder',linspecer(length(OUT.K))*.5) ;
                plot([1;1]*real(OUT.K),get(gca,'ylim').','linewidth',1)
                set(gca,'colororderindex',1) ;
                plot([1;1]*real(OUT.K+OUT.dK),get(gca,'ylim').','-.','linewidth',1)
                plot([1;1]*real(OUT.K-OUT.dK),get(gca,'ylim').','-.','linewidth',1)




%%
K = OUT.K ;
U = OUT.U ;
ESTER = OUT.ESTER ;
f = K/dt/2/pi ;
[~,indSort] = sort(real(f)) ;


clf ;
    freq = (0:nT-1)/nT*Fe ;
    plot(freq,meanUx,'linewidth',2)
    plot(freq,meanUy,'linewidth',2)
    plot(freq,meanUz,'linewidth',2)
    plot(freq,meanVel,'-.k','linewidth',1)
    axis tight
    set(gca,'xlim',[10 15000])
    set(gca,'xscale','log','yscale','log')
    grid on
    %legend({'$|\dot u_x|$','$|\dot u_y|$','$|\dot u_z|$','$||\dot \mathbf{u}||$'},'edgecolor','k')
    xlabel('frequency (Hz)','units','normalized','position',[1 0],'verticalalignment','bottom','horizontalalignment','right')
    plot([real(f);real(f)],repmat(get(gca,'ylim')',[1 length(f)]),':k','linewidth',1)
    plot(real(f),imag(f)./real(f)/10,'.k','markersize',25) ;
    
% SAVE THE RESULTS OF MODAL ANALYSIS
    [file,path] = uiputfile('_HDMA.mat','SAVE THE MODAL ANALYSIS RESULTS') ;
    if file~=0
        % Sort frequencies
            [~,indSort] = sort(real(f),'ascend') ;
            modalFreq = f(indSort) ;
            modalU = U(:,:,indSort) ;
        % Keep positive only
            indPos = find(real(modalFreq)>0 & imag(modalFreq)>0) ;
            modalFreq = modalFreq(indPos) ;
            modalU = modalU(:,:,indPos) ;
        save([path,file],'X0','Y0','Z0','modalFreq','modalU','indP') ;
    end
    
    

%% MEAN VELOCITY FIGURE

clf ;
    freq = (0:nT-1)/nT*Fe ;
    plot(freq,meanUx,'linewidth',2)
    plot(freq,meanUy,'linewidth',2)
    plot(freq,meanUz,'linewidth',2)
    plot(freq,meanVel,'-.k','linewidth',1)
    axis tight
    set(gca,'xlim',[10 1000])
    set(gca,'xscale','log','yscale','log')
    grid on
    legend({'$|\dot u_x|$','$|\dot u_y|$','$|\dot u_z|$','$||\dot \mathbf{u}||$'},'edgecolor','k')
    xlabel('frequency (Hz)','units','normalized','position',[1 0],'verticalalignment','bottom','horizontalalignment','right')
    ylim = get(gca,'ylim') ;
    ff = real(sort(f)) ;
    ff = ff(ff>0) ;
    linescolor = [1 1 1]*0.3 ;
    plLines = plot([real(f);real(f)],repmat(get(gca,'ylim')',[1 length(f)]),'-','color',linescolor,'linewidth',1) ;
    uistack(plLines,'bottom') ;
    txt = [] ;
    for i = 1:length(ff)
        txt(i) = text(real(ff(i)),ylim(1),['\bf \,',char(i+96),'.'],'horizontalalignment','left','verticalalignment','bottom','color',linescolor) ;
    end
    %set(txt,'edgecolor','k')
    
    
%% LIEP INITIALIZATION

    cd 'C:\Users\pierre.margerit\Documents\MATLAB\THESE\+DATA\Mesures\Robot Laser\Poutre_S\ARTICLE'
    %EXP = load('Curved2_HDMA') ;
    %EXP = load('Curv100Pts_HDMA') ;
    EXP = load('300Pts_newESPRIT_HDMA') ;
    FEM = load('Curved_ElemMat') ;

    % EXTRACT INFOS
        EXP.X0 = EXP.X0(EXP.indP) ;
        EXP.Y0 = EXP.Y0(EXP.indP) ;  
        EXP.Z0 = EXP.Z0(EXP.indP) ;
        FEM.noCL_NODES = FEM.NODES(FEM.indNoCL_U,:) ;
        EXP.nF = length(EXP.modalFreq) ;
        EXP.nPts = length(EXP.X0) ;
        FEM.nPts = size(FEM.noCL_NODES,1) ;
        
    % ROTATION DES POUTRES POUR SE METTRE DANS LE MEME REPERE
        scale = 1 ;
        % Translation et Rotation
        FEM.AbsNODES = FEM.noCL_NODES-repmat(FEM.NODES(1,:),[FEM.nPts 1]) ;
        EXP.AbsX0 = -(EXP.X0-EXP.X0(end))/1000*scale ;
        EXP.AbsY0 = (EXP.Z0-EXP.Z0(end))/1000*scale ;
        EXP.AbsZ0 = (EXP.Y0-EXP.Y0(end))/1000*scale ;
        EXP.AbsNODES = [EXP.AbsX0,EXP.AbsY0,EXP.AbsZ0] ;

        EXP.U = [-EXP.modalU(:,1,:) EXP.modalU(:,3,:) EXP.modalU(:,2,:)] ;

        cla ;
            plot3(FEM.AbsNODES(:,1),FEM.AbsNODES(:,2),FEM.AbsNODES(:,3),'k')
            plot3(EXP.AbsX0,EXP.AbsY0,EXP.AbsZ0,'.r')
            axis equal 
            xlabel 'X'
            ylabel 'Y'
            zlabel 'Z'
    
    % MATRICE D'OBSERVATION C (POINT LE PLUS PROCHE)
            FEM.l = sqrt(sum(diff([0 0 0 ; FEM.AbsNODES],1,1).^2,2)) ;
            L = max(FEM.l)*2 ;
        % Interpolation d'un DDL par la paramétrage des courbes
            C = sparse(EXP.nPts,FEM.nPts) ;
            pMissed = [] ;
            for pt = 1:EXP.nPts
                dist = sqrt(sum((repmat(EXP.AbsNODES(pt,:),[FEM.nPts 1])-FEM.AbsNODES).^2,2)) ;
                nearest = find(dist==min(dist),1) ;
                if dist(nearest)>0.51*L ; pMissed(end+1) = pt ; continue ; end
                C(pt,nearest) = 1 ;
            end
            pMissed 
        % Assemblage de tous les DDLs
            if 1 % Add the shift (measurements are made on the top surface)
                h = 0.004 ;
                C = [C*1 C*0 C*0 C*h C*0 C*0 ; ...
                     C*0 -C*h C*1 C*0 C*0 C*0 ; ...
                     C*0 C*0 C*0 C*0 C*1 C*0 ; ...
                    ] ;
            else
                C = [C*1 C*0 C*0 C*0 C*0 C*0 ; ...
                     C*0 C*0 C*1 C*0 C*0 C*0 ; ...
                     C*0 C*0 C*0 C*0 C*1 C*0 ; ...
                    ] ;
            end
            C(:,[1:6:end,2:6:end,3:6:end,4:6:end,5:6:end,6:6:end]) = C  ;
            
            clf ; plotMatrix(log10(abs(C))) ; colormap(flip(gray))
            
            
%% LIEP POUTRE ISOTROPE
            
    cullModes = [] ;
    
    indModes = setdiff(1:EXP.nF,cullModes) ;
    
    vararg = {'solver','eigs'...
                'stepRatio',1 ...
            	'Nfa',round(1.3*length(indModes)) ...
                'weightPowF',-.5 ...
                'maxIt',1000 ...
                'absTolF',1e-2 ...
                'relTolF',1e-3 ...
                'tolE',1e3 ...
                'monitoring',true ...
                'video',false ...
                'corrMin',0.6 ...
                'paramNames',{'E','G'} ...
                ...'K0',FEM.K_e*Eei ...
                } ;
    
    Ei = 3.5e9 ;
    nui = .3 ;
    Gi = Ei/2/(1+nui) ;
    rho = 1190 ;
    e0 = [Ei,Gi] ;
    
    f = (EXP.modalFreq) ;
    w2 = (2*pi*f).^2 ;
    U = reshape(EXP.U,[],EXP.nF) ;
    K = {FEM.K_b1+FEM.K_b2+FEM.K_e,FEM.K_s1+FEM.K_s2+FEM.K_t} ;
    M = FEM.M*rho ;
    
    clc
    close all
    indModes = setdiff(1:EXP.nF,cullModes) ;
    [e,out] = LIEP(w2(indModes),U(:,indModes),C,K,M,e0,vararg{:}) 
    nu = e(1)/2/e(2) - 1
    
    
%% FIGURE MACS

    margin = 10 ;
    titleHeight = 70 ;
    
    MAC = abs(out.MAC(:,:,end)) ;
    close all ;
    fig = figure() ;
    fig.Position(3:4) = size(MAC')*50 ;
    axes('outerposition',[0 0 1 1]) ;
    plotMatrix(MAC) ;
    colormap(flip(gray)) ;
    txt = [] ;
    for i = 1:size(MAC,1)
        for j = 1:size(MAC,2)
            txt(end+1) = text(j,i,1,num2str(MAC(i,j)*100,'%.1f')) ;
            if MAC(i,j)>0.5
                set(txt(end),'Color',[1 1 1]*.99) ;
            end
        end
    end
    set(gca,'xtick',[1:size(MAC,2)]','xticklabel',strcat('\bf ',num2str([1:size(MAC,2)]'))) ;
    set(gca,'ytick',[1:size(MAC,1)],'yticklabel',strcat('\bf{',char(96+(1:size(MAC,1)))','}~')) ;
    set(txt,'fontsize',15) ;
    xlabel('\bf Numerical Modes')
    ylabel('\bf Experimental Modes')
    set(gca,'sortmethod','childorder')
    
    
%% PLOT RESULTS OF MATCHING

        nIt = out.it ;
        it = nIt ;
        amp = 0.1 ;
        
        nA = size(out.MAC,1) ;
        nAx = ceil(sqrt(nA)) ;
        nAy = ceil(nA/nAx) ;
        X0 = EXP.AbsX0 ;
        Y0 = EXP.AbsY0 ;
        Z0 = EXP.AbsZ0 ;
        Lx = max(X0(:))-min(X0(:)) ; 
        Ly = max(Y0(:))-min(Y0(:)) ;
        Lz = max(Z0(:))-min(Z0(:)) ;
        NORM = norm([Lx Ly Lz]) ;
    
    close all
    fig = figure ;
        AXES = [] ;
        S0 = [] ;
        Sf = [] ;
        Sfa = [] ;
        Titles = [] ;
        for ay = 1:nAy
            for ax = 1:nAx
                ff = (ay-1)*nAx+ax ;
                if ff>nA ; break ; end 
                AXES(ff) = mysubplot(nAy,nAx,ff) ;
                    % First infos
                        ttlStr = {} ;
                        fnum = out.matchIndF(ff,it) ;
                        %ttlStr{end+1} = ['Exp : ',complex2str(f(ff)),' Hz'] ;
                    % Neutral Axis
                        S0(ff) = plot3(X0,Y0,Z0,'--k','linewidth',.5) ;
                    % Experimental shape
                        Uexp = EXP.U(:,:,ff) ;
                        normU = sqrt(sum(abs(Uexp).^2,2)) ;
                        [maxExpU,ptMax] = max(normU) ;
                        [~,compMax] = max(abs(Uexp(ptMax,:))) ;
                        phase = exp(1i*angle(Uexp(ptMax,compMax))) ;
                        AmpExp = amp*NORM/maxExpU/phase ;
                        Sf(ff) = plot3(X0+real(AmpExp*Uexp(:,1)),...
                                        Y0+real(AmpExp*Uexp(:,2)),...
                                        Z0+real(AmpExp*Uexp(:,3)),...
                                        'ok','markersize',5,'linewidth',.5) ;
                    % Model-given shape
                        if fnum~=0  % matched freq.
                            %ttlStr{end+1} = ['Num : ',complex2str(sqrt(out.fa(fnum,it))/2/pi),' Hz'] ;
                            ttlStr{end+1} = ['MAC : ',num2str(abs(out.MAC(ff,fnum,it))*100,'%.2f'),' %'] ;
                            % On the measurement points
                                Ua = reshape(out.Ua(:,fnum,it),[],3) ;
                                Ua = Ua/out.MAC(ff,fnum,it)/phase ;
                                maxNumU = max(sqrt(sum(abs(Ua).^2,2))) ;
                                AmpNum = amp*NORM/maxNumU ;
                                Sfa(ff) = plot3(X0+real(AmpNum*Ua(:,1)),...
                                                Y0+real(AmpNum*Ua(:,2)),...
                                                Z0+real(AmpNum*Ua(:,3)),...
                                                 '-r','linewidth',2) ;
                            % On the model mesh
                                if 0
                                    Va = reshape(out.Va(:,fnum),6,[]).' ;
                                    Va = Va(:,1:2:5)/out.MAC(ff,fnum,it) ;
                                    maxNumV = max(sqrt(sum(real(Va).^2,2))) ;
                                    AmpV = amp*NORM/maxNumV ;
                                    plot3(FEM.AbsNODES(:,1)+AmpV*real(Va(:,1)),...
                                            FEM.AbsNODES(:,2)+AmpV*real(Va(:,2)),...
                                            FEM.AbsNODES(:,3)+AmpV*real(Va(:,3)),...
                                            '-b','linewidth',1) ;
                                end
                        end
                        % Frequency titles
                            Titles(ff) = title(ttlStr,...
                                                ...'position',[0.5 1 0],...
                                                ...'horizontalalignment','center',...
                                                ...'verticalalignment','top',...
                                                'interpreter','none',...
                                                'fontweight','normal',...
                                                'fontsize',12) ;
            end
        end
        axis(AXES,'equal') ;
        axis(AXES,'tight') ;
        axis(AXES,'off') ;
        set(AXES,'view',[40 22])
        set(AXES,'xtick',[],'ytick',[],'ztick',[])
        global hlink
        hlink = linkprop(AXES,'view') ;
        
        
%% SAVE ALL AXES AS INDIVIDUAL FIGURES

    fig = figure ;
    set(AXES,'units','pixels') ;
    figSize = get(AXES(1),'position') ;
    set(AXES,'units','normalized') ;
    fig.Position(3:4) = figSize(3:4)*1.5 ;
    for a = 1:length(AXES)
        clf ;
        ax = copyobj(AXES(a),fig) ;
        ax.Position = [0 0 1 1] ;
        ax.XTick = [] ;
        ax.YTick = [] ;
        ax.Title = [] ;
        ax.CameraViewAngle = ax.CameraViewAngle;%/1.9 ;
        if a == 1
            disp(['RESIZE AND ROTATE THE VIEW, THEN PRESS ENTER']) ;
            pause() ;
            if ~isvalid(fig) ; return ; end
            camViewAngle = ax.CameraViewAngle ;
            camPos = ax.CameraPosition ;
            camTarget = ax.CameraTarget ;
        end
        ax.CameraViewAngle = camViewAngle ;
        ax.CameraPosition = camPos ;
        ax.CameraTarget = camTarget ;
        drawnow ;
        recfig(['Figures/CompModes/Mode',num2str(a)]) ;
    end
    
    



