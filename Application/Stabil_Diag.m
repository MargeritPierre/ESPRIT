%% ===========================================================
% STABILIZATION DIAGRAM
% ============================================================

% -----------------------------------------------------
    % DATA TO PROCESS
        U = [AvgdX ; AvgdY ; AvgdZ] ; % Total Measurement Data [x,t]
        indT = indT ; % Time indices used for ESPRIT
        OUT = OUT ; % ESPRIT results
        Fe = Fe ; %Sampling Frequency
        XLIM = [10 8000] ;
        
    % Plot a specific Criterion
        process =   { ...
                      ... 'sorting' ... % Sort by MAC values
                      ... 'stability' ... % Stability of the pole between successive orders
                      ... 'complexity' ... % Measure out-of phase behavior
                       'histogram' ... % Make an histogram
                    };
        vizualize = true ; % some tools with the mouse
    % Criterion parameters
        % Common
            %colors = jet(1000)*.85 ;
            colors = bsxfun(@plus,bsxfun(@times,gray(1000),[0 1 1]),[1 0 0]) ;
        % Sorting
            minMAC = 0 ; 
        % Complexity
            intervalComplex = [0 1] ; % for color axis
        % Stability
            paramStab =  'frequency' ... % Complex frequency
                        ... 'mode' ... % Mode
                        ;
            maxRelDiff = 0.005 ; % Maximum relative variation allowed between two successive orders
            intervalStab = [1e-4 1] ;
        % Histogram
            maxRelSTD = 1e-3 ; % Maximum relative standard deviation, used to build the edges.
% -----------------------------------------------------
    
            
    % PROCESSING        
    % Signal Order infos
        R0 = OUT.R0;
        minOrder = min(R0) ;
        maxOrder = max(R0) ;
        
    % FORMAT MEASUREMENTS SPECTRUM
        fftU = sum(abs(fft(U,[],2)).^2,1) ;
        freqU = (0:length(fftU)-1)/length(fftU)*Fe ;
        modfftU = log10(fftU) ;
        modfftU = (maxOrder-minOrder)*(modfftU-min(modfftU(:)))/(max(modfftU(:))-min(modfftU(:)))+minOrder ;

    % FORMAT FREQUENCIES
        Kstab = OUT.Kstab*Fe/2/pi ;
        Ustab = OUT.Ustab ;
    % Cull negative frequencies
        indPos = real(Kstab)>0 ;
        Kstab(~indPos) = NaN*(1+1i) ;
        for n = 1:size(Ustab,1) ;
            uu = Ustab(n,:,:) ;
            uu(~indPos) = NaN*(1+1i) ;
            Ustab(n,:,:) = uu ;
        end
    % Re-sort frequencies
        for rr=1:size(Kstab,1) ;
            [Kstab(rr,:),indsort] = sort(Kstab(rr,:),'ascend') ;
            Ustab(:,rr,:) = Ustab(:,rr,indsort) ;
        end
        allNaN = ~any(~isnan(Kstab),1) ;
        Kstab = Kstab(:,~allNaN) ;
        Ustab = Ustab(:,:,~allNaN) ;
        
    % INIT FIGURE
        clf('reset')
        ax = gca ;
            ptsK = plot(abs(real(Kstab)),R0,'.k','markersize',8) ;
            %ptsK = scatter(abs(real(Kstab(:))),repmat(R0(:),[size(Kstab,2) 1]),10,'k') ;
            spectrum = plot(freqU,modfftU,'-','linewidth',2,'color',[1 1 1]*.5) ;
            ax.XScale = 'log' ;
            if isempty(XLIM)
                ax.XLim(2) = Fe/2 ;
            else
                ax.XLim = XLIM ;
            end
            mark = plot(NaN,NaN,'ob') ;
            drawnow ;
            grid on
        zoom('xon') ;
    
    % PROCESS THE RESULTS
        for i = 1:length(process)
            switch process{i}
                case 'sorting' % MAC VALUES
                    % Init plot
                        branches = plot(NaN,NaN,'-r','linewidth',1) ;
                    % Init Data
                        KstabSort = NaN*(1+1i)*ones(size(Kstab)) ;
                        KstabSort(end,:) = Kstab(end,:) ;
                        UstabSort = NaN*(1+1i)*ones(size(Ustab)) ; 
                        UstabSort(:,end,:) = Ustab(:,end,:) ;
                    % Build Branches
                        for or = length(R0)-1:-1:1 ;
                            % Get the modes
                                U1 = squeeze(UstabSort(:,or+1,:)) ;
                                U2 = squeeze(Ustab(:,or,:)) ;
                            % Compute the MAC matrix
                                mac = abs(U2'*U1)./sqrt(sum(abs(U2).^2,1)'*sum(abs(U1).^2,1)) ;
                            % Match the modes
                                oldInd = [] ; newInd = [] ;
                                for rr = 1:size(mac,1)
                                    % Find the max
                                        if ~any(~isnan(mac(:))); break; end
                                        maxMAC = max(mac(:)) ;
                                        if maxMAC<minMAC ; break ; end
                                        [rmax,cmax] = find(mac==maxMAC) ;
                                        oldInd(end+1) = cmax(1) ;
                                        newInd(end+1) = rmax(1) ;
                                    % Set the corresponding values to NaN
                                        mac(:,cmax(1)) = nan ;
                                        mac(rmax(1),:) = nan ;
                                end
                                % Save
                                    KstabSort(or,oldInd) = Kstab(or,newInd) ;
                                    UstabSort(:,or,oldInd) = U2(:,newInd) ;
                                % Plot branches
                                    branches.XData = abs(real(KstabSort(:))) ;
                                    branches.YData = repmat(R0(:),[size(Kstab,2) 1]) ;
                                    drawnow ;
                        end
                    % Re-Plot individual branches
                        delete(branches) ;
                        ax.ColorOrder = [1 1 1]*0+[1 0 0]*1 ;  jet(size(KstabSort,2)) ;
                        branches = plot(abs(real(KstabSort)),R0(:),'-','linewidth',1) ;
                    % Set the points on top
                        pts2 = copyobj(ptsK,gca) ;
                        delete(ptsK) ;
                        ptsK = pts2 ;
                        clear('pts2') ;
                        drawnow ;
                case 'stability'
                    if ~ismember(process,'sorting'); warning(['Stability cannnot be measured without sorting']) ; end
                    % Compute relative variations
                        switch paramStab
                            case 'frequency'
                                Var = [diff(KstabSort,1) ; zeros(1,size(KstabSort,2))] ;
                                relVar = Var./KstabSort ;
                                outOfVarMax = abs(relVar)>maxRelDiff ;
                                critStab = [abs(diff(abs(KstabSort),1,1))./(abs(KstabSort(1:end-1,:))+abs(KstabSort(2:end,:))) ; NaN*ones(1,size(KstabSort,2))] ;
                                %critStab = [abs(diff(abs(KstabSort),1,1)) ; NaN*ones(1,size(KstabSort,2))] ;
                            case 'mode'
                        end
                    % Segregate frequencies
                        KstabVar = KstabSort ; 
                        KstabVar(outOfVarMax) = NaN*(1+1i) ;
                    % Plot
                        if 0 % Newbranches
                            delete(branches)
                            %plot(abs(real(KstabVar)),R0(:),'-k','linewidth',2)
                        else % Criterion Plot
                            set(branches,'color',[1 1 1]*.7) ;
                            delete(branches)
                            delete(ptsK) ;
                            ptsK = scatter(abs(real(KstabSort(:))),repmat(R0(:),[size(KstabSort,2) 1]),10,log10(critStab(:)),'filled') ;
                            colorbar
                            colormap(colors)
                            caxis(log10(intervalStab)) ;
                        end
                case 'complexity' % MODE COMPLEXITY
                    complex = zeros(size(Kstab))*NaN ;
                    for ii = 1:size(Kstab,1)
                        for jj = 1:size(Kstab,2)
                            if isnan(Kstab(ii,jj)); continue; end
                            uu = Ustab(:,ii,jj) ;
                            ss = [real(uu(:)) imag(uu(:))] ;
                            s = sqrt(eig(ss'*ss,'vector')) ;
                            complex(ii,jj) = (min(s))/max(s) ;
                        end
                    end
                    delete(ptsK) ;
                    ptsK = scatter(abs(real(Kstab(:))),repmat(R0(:),[size(Kstab,2) 1]),10,complex(:),'filled') ;
                    colorbar
                    colormap(colors)
                    caxis(intervalComplex)
                case 'histogram'
                    % Prepare the axes
                        axHist = axes('position',ax.Position);
                        axHist.XLim = XLIM ;
                        axHist.XTick = [] ;
                        axHist.XScale = 'log' ; 
                        axHist.YAxisLocation = 'right' ;
                    % Make the histogram
                        edges = 10.^(log10(XLIM(1)):maxRelSTD:log10(XLIM(2))) ;
                        %edges = logspace(log10(XLIM(1)),log10(XLIM(2)),hist_count+1) ;
                        [numFreq,edges,indices] = histcounts(real(Kstab(:)),edges) ;
                        if 1
                            histHandle = histogram(real(Kstab(:)),edges) ;
                            histHandle.EdgeColor = 'r' ;
                            histHandle.FaceColor = 'r' ;
                        else
                            histHandle = plot((edges(1:end-1)+edges(2:end))/2,numFreq,'r','linewidth',1)
                        end
                        set(ptsK,'markersize',4) ;
                    % Link the axes
                        uistack(ax,'top') ;
                        linkaxes([axHist ax],'x') ;
            end
        end
        
        
    % INTERACTIVE VIZUALIZATION
        if vizualize
            % Init plot
                axShape = axes() ;
                    axShape.Position = [ax.Position(1:2) .2 .2] ;
                    box on
                    shapeR = plot(axShape,real(Ustab(:,end,end))*0+1,'b') ;
                    shapeI = plot(axShape,imag(Ustab(:,end,end))*0,':b') ;
                    axShape.YLim = [-1 1] ;
                    axShape.YTick = [] ;
                    axShape.XTick = [] ;
                axPolar = axes() ;
                    axPolar.Position = [axShape.Position(1:2)+axShape.Position(3:4).*[1 0] .2 .2] ;
                    box on
                    axis equal
                    axPolar.XLim = [-1 1] ;
                    axPolar.YLim = [-1 1] ;
                    axPolar.YTick = [] ;
                    axPolar.XTick = [] ;
                    polarShape = plot(axPolar,NaN,NaN,'+b','linewidth',.5,'markersize',5) ;
            % Loop
                while true
                    % Get the point
                        cp = ax.CurrentPoint(1,1:2) ;
                        ord = closest(R0,round(cp(2))) ;
                        kk = closest(abs(real(Kstab(ord,:))),cp(1)) ;
                        disp(char(10)) ; disp(['order: ',num2str(ord)]) ; disp(['num: ',num2str(kk)]) 
                    % Get the mode
                        uu = Ustab(:,ord,kk) ;
                        %uu = uu/sqrt(mean((uu./abs(uu)).^2)) ;
                        uu = uu/max(abs(uu(:))) ;
                    % Update the figure
                        mark.XData = abs(real(Kstab(ord,kk))) ;
                        mark.YData = R0(ord) ;
                        shapeR.YData = real(uu) ;
                        shapeI.YData = imag(uu) ;
                        polarShape.XData = real(uu) ;
                        polarShape.YData = imag(uu) ;
                        drawnow ;
                end
        end
    
    
%% PROCESS THE MODES

% S_BEAM
    ORD_NUM=[100,1 ; 100,2 ; 100,4 ; 100,6 ; 100,8 ; 100,10 ; 100,15 ; 37,14 ; 34,15 ; 30,14 ; 37,19 ; 35,19 ; 46,28 ; 45,29 ; 45,30 ; 54,39 ; 54,40 ; 67,55 ; 64,53 ; 55,45 ; 64,55 ; 67,59 ; 69,62 ; 68,62 ; 73,68 ; 70,66 ; 74,71 ; 100,98 ; 100,99] ;
    NPtsY = 3 ;
    
% GET THE INFOS
    ORD = ORD_NUM(:,1) ;
    NUM = ORD_NUM(:,2) ;
    fOK = diag(Kstab(ORD,NUM)) ;
    
% SIGNAL MODEL
    t = (0:size(U(:,indT),2)-1)/Fe ;
    Vand = [exp(2i*pi*t'*fOK(:).') exp(-2i*pi*t'*conj(fOK(:)).')] ;
    A = Vand\U(:,indT).' ;
    SignalModel = (Vand*A).' ;
    fftSM = abs(fft(SignalModel,[],2)).^2 ;
    fftU = abs(fft(U(:,indT),[],2)).^2 ;
    freqU = (0:length(fftU)-1)/length(fftU)*Fe ;

% RESHAPE
    fftU = squeeze(sum(reshape(fftU,[],NPtsY,length(indT)),1)).' ;
    fftSM = squeeze(sum(reshape(fftSM,[],NPtsY,length(indT)),1)).' ;
    
% INIT FIGURE
    clf('reset')
    ax = gca ;
        plot(freqU,fftU,'-','linewidth',1) ;
        set(gca,'colororderindex',1)
        plot(freqU,fftSM,':','linewidth',2) ;
        ax.XScale = 'log' ;
        ax.YScale = 'log' ;
        axis tight
        if isempty(XLIM)
            ax.XLim(2) = Fe/2 ;
        else
            ax.XLim = XLIM ;
        end
        grid on
% Plot frequency info
    %plot([1;1]*real(fOK(:)).',get(gca,'ylim'),'-.k','linewidth',1)
    

%% Frequency and damping
    clf
    plot(real(fOK),imag(fOK)./real(fOK),'o-')
    set(gca,'xscale','log','yscale','log')
    grid on
    set(gca,'ylim',[0.01 0.05])
    
    
%% SAVE THE FIGURE AS PDF
    NAME = 'Poutre_S_StabDiag_Histogram'
    NptsMax = 2000 ;
    
    % Compress the spectrum data
    xspec = spectrum.XData ;
    yspec = spectrum.YData ;
    XLIM = get(gca,'xlim') ;
    xcomp = unique(logspace(log10(min(XLIM)),log10(max(XLIM)),NptsMax)) ;
    ycomp = interp1(xspec,yspec,xcomp) ;
    spectrum.XData = xcomp ;
    spectrum.YData = ycomp ;
    
    % Save
    axs = findall(gcf,'type','axes') ;
    pos = get(axs,'outerposition') ;
    set(axs,'position',[0 0 1 1]) ; axis(axs,'off')
    drawnow ;
    recfig(NAME)
    for aa = 1:length(axs)
        set(axs(aa),'outerposition',pos{aa}) ; axis(axs,'on')
    end
    
    % Reset compression
    spectrum.XData = xspec ;
    spectrum.YData = yspec ;
    drawnow ;



    
    
    
