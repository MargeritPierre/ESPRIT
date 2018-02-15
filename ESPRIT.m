function [K,U,ESTER] = ESPRIT(Signal,varargin)
%
% ESPRIT Extract wavevectors from a signal
%
%   varargout = ESPRIT(Signal,varargin)
%   [K,U,ESTER] = ESPRIT(Signal,'Name','Value',...)
%   [K,U,ESTER] = ESPRIT(Signal,DIMS_K,'Name','Value',...)
%
%   varargin parameters (and default values) : 
%       CHOICE : 'auto' or 'manual' ('auto')
%       ESTER_THRS : ester criterion threshold (1)
%       DIMS_K : dimensions of the wavevectors (ndims(Signal))  
%       FUNC : type of function searched 'exp' or 'cos' (repmat('exp',[length(DIMS_K) 1]))
%       R0 : signal order candidates (1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)))
%       FIT : 'LS' or 'TLS' (TLS)s
%       DECIM : decimation factors (ones(1,ndims(Signal)))
%           - along DIMS_K : /!\ NYQUIST
%           - along DIMS_P : reduce the size of Css. However, U is estimated at all points.
%       SHIFTS : for multiresolution (eye(length(DIMS_K)))
%       DEBUG : bool to prompt procedure state (false)
%
%   varargout : 
%       K : extracted complex wavevectors (size : [length(DIMS_K) R])
%       U : complex amplitudes (modes) (size : [size(Signal(DIMS_P)) R])
%       ESTER : ester criterion for each shift (size : [size(SHIFTS,1) length(R0)])
%
%
% This ESPRIT algorithm implements :
%   - Standard ('exp') and Paired ('cos') ESPRIT (FUNC)
%   - ESTER criterion (ESTER_THRS)
%   - Multidimensionnal (DIMS_K)
%   - LS or TLS linear regression (FIT)
%   - Multiple Invariance or Multiple Resolution (SHIFTS)
%   - Decimate ESPRIT (DECIM)
%    
% To go further on developments :
%   (1) use of high-order cummulents
%   (2) MR-ESPRIT : solving ambiguity when Nyquist is violated
%   (3) GUI for manual choice
%   (4) Cramer-Rao bound, uncertainty estimation
%   (5) replace kron() by repelem() in indHbH_d build for speed (ver>=Matlab2015)
%   (6) choice between EIG, EIGS, SVD or SDS for computing W for speed


% INITIALIZATION

    % Input Initialization
        parseInputs ; 
        varargin ;
        CHOICE ;
        ESTER_THRS ; 
        DIMS_K ; 
        FUNC ; 
        R0 ; 
        FIT ; 
        DECIM ; 
        SHIFTS ; 
        DEBUG ;

   % MANUAL INPUT MODIF. (for debugging)
        %SHIFTS = [1 1 1 ; 0 -1 1 ; 0 1 1] ;
        %SHIFTS = [SHIFTS ; 2*SHIFTS] ;
        %SHIFTS = [1 2 3]' ;
        %DECIM = [1 1 2] ;
        %FUNC = ['exp' ; 'cos' ; 'exp'] ;
        %FUNC = repmat('cos',[length(DIMS_K) 1]) ;
        
   % DEBUGGING ?     
        if(DEBUG)
            display(char(10)) ;
            display('----------- ESPRIT ND -----------') ;
            display('   Initialization : ') ; lastTime = tic ;
        end
        
    % Necessary informations
        SIZE = size(Signal) ; % Size of the signal
        NDIMS = ndims(Signal) ; % Dimension of the signal
        DIMS_P = setdiff(1:NDIMS,DIMS_K) ; % Dimensions of the isosurface
        Lk = SIZE(DIMS_K) ; % Number of points in the dimensions of the exponentials
        Lp = SIZE(DIMS_P) ; % Number of points in the dimensions of the isophase surface
    
    % Reshaping the signal
        Signal = permute(Signal,[DIMS_P DIMS_K]) ; % Sort dimensions
        Signal = reshape(Signal,[prod(Lp) prod(Lk)]) ;
        
    % Output Initialization
        K = zeros([max(R0) length(DIMS_K)]) ;
        U = zeros([max(R0) length(DIMS_P)+1]) ;
        ESTER = zeros([size(SHIFTS,1) length(R0)]) ;
        
        
        
        
% HANKEL MATRIX BUILDING
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Hankel Matrix Construction : ') ; lastTime = tic ; end
    
    % Sub-Hankel matrices shapes (Kk rows, Mk columns)
        DECIM_K = DECIM(DIMS_K) ;
        Kk = zeros([1 length(DIMS_K)]) ;
        Mk = zeros([1 length(DIMS_K)]) ;
        isCOS = false([1 length(DIMS_K)]) ; % IS COSINUS SEARCHED ?
        for d = 1:length(DIMS_K)
            switch lower(FUNC(d,:))
                case 'exp'
                    Kk(d) = floor((Lk(d)+DECIM_K(d))./(1+DECIM_K(d))) ; % max size of (signal+noise) subspace
                    Mk(d) = Lk(d) - DECIM_K(d).*(Kk(d)-1) ;
                case 'cos' % /!\ decimation not included !
                    Kk(d) = floor((Lk(d)+1)/3) ;
                    Mk(d) = Lk(d)-2*Kk(d)+1 ;
                    isCOS(d) = true ;
            end
        end
        
    % Adjust signal order candidates
        R0 = R0(R0<=prod(Kk)-1-max(sum(abs(SHIFTS),2))) ;
        
    % Sub-Hankel matrices indices
        subIndH = cell(length(DIMS_K),1) ;
        for d = 1:length(DIMS_K)
            switch lower(FUNC(d,:))
                case 'exp'
                    subIndH{d,1} = (1:DECIM_K(d):DECIM_K(d)*(Kk(d)-1)+1)'*ones(1,Mk(d)) + ones(Kk(d),1)*(0:Mk(d)-1) ;
                case 'cos' % /!\ decimation not included !
                    subIndH{d,1} = flipud(hankel(0:Kk(d)-1,Kk(d)-1:Lk(d)-Kk(d)-1)) + 1 ;
                    subIndH{d,2} = hankel(Kk(d):2*Kk(d)-1,2*Kk(d)-1:Lk(d)-1) + 1 ;
            end
        end
        
    % Hankel-block-Hankel matrix indices
        
        % Indices for each dimension (used in shift invariances and Vandermonde matrix)
            indHbH_d = cell(length(DIMS_K),1) ;
            D = 1:length(DIMS_K) ;
            for d = 1:length(D)
                switch lower(FUNC(d,:))
                    case 'exp'
                        indHbH_d{d,1} = repmat(subIndH{d,1},[prod(Kk(D<d)), prod(Mk(D<d))]) ;
                        indHbH_d{d,1} = kron(indHbH_d{d,1},ones([prod(Kk(D>d)), prod(Mk(D>d))])) ;
                    case 'cos'
                        for i = 1:2
                            indHbH_d{d,i} = repmat(subIndH{d,i},[prod(Kk(D<d)), prod(Mk(D<d))]) ;
                            indHbH_d{d,i} = kron(indHbH_d{d,i},ones([prod(Kk(D>d)), prod(Mk(D>d))])) ;
                        end
                end
            end
            
        % BUILD THE MATRIX
            H = buildHbH ;
        
        
        
% AUTOCOVARIANCE MATRIX
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Computing of Css : ') ; lastTime = tic ; end

        Css = H*H'/(prod(Mk)*prod(Lp)) ;
        
        
        
% EIGENVALUE DECOMPOSITION        
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Eigendecomposition of Css : ') ; lastTime = tic ; end
    % EIG
        [W,D] = eig(Css,'vector') ;
        [~,Ind] = sort(D,'descend') ;
        W = W(:,Ind) ;
        
        
% ESTER CRITERION
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   ESTER Criterion : ') ; lastTime = tic ; end
    % Compute all errors
        ESTER = zeros(size(SHIFTS,1),length(R0)) ;
        for r = 1:length(R0)
            for s = 1:size(SHIFTS,1)
                ESTER(s,r) = ester(r,s) ;
            end
        end
    % Compute the mean
        %CRIT = prod(ESTER,1) ;
        CRIT = max(ESTER,[],1) ;
        %CRIT = mean(ESTER,1) ;
    % CHOOSE THE ORDER
        R = find(CRIT>=max(CRIT)*ESTER_THRS,1,'last') ;
    
% STABILIZATION DIAGRAM
    if STABILDIAG
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Stabilization Diagram : ') ; lastTime = tic ; end
        stabilizationDiagram() ;
    end
        
    
% EXTRACT THE WAVEVECTORS
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Pole Estimation : ') ; lastTime = tic ; end
    if (nargout>=1) ; extractPoles(R) ; end
        
    
% ESTIMATE AMPLITUDES
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Amplitude Estimation : ') ; lastTime = tic ; end
    if (nargout>=2) ; computeU ; end
    

% END OF THE PROCEDURE
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('---------------------------------') ; end


    
    
% ===================================================================================================================    
% FUNCTIONS FOR ESPRIT
% ===================================================================================================================

% HANKEL-BLOCK-HANKEL MATRIX
    function H = buildHbH
        % Decimation in the isosurface points
            indP = 1:size(Signal,1) ;
            DECIM_P = DECIM(DIMS_P) ;
            if any(DECIM_P~=1)
                indP = (1:DECIM_P(1):Lp(1))' ;
                D = 1:length(DIMS_P) ;
                for d = 2:length(DECIM_P)
                    indP_d = (1:DECIM_P(d):Lp(d))' ;
                    indP = repmat(indP,[length(indP_d) 1]) + kron((indP_d-1)*prod(Lp(D<d)),ones(length(indP),1)) ;
                end
            end
        % HbH INDICES
            if size(indHbH_d,2)==1 % ONLY EXPONENTIALS SEARCHED
                n_indHbH = 1 ;
                indHbH{1} = indHbH_d{1,1} ;
                D = 1:length(DIMS_K) ;
                for d = 2:length(DIMS_K)
                    indHbH{1} = indHbH{1} + (indHbH_d{d}-1)*prod(Lk(D<d)) ;
                end
            else % SOME COSINUS SEARCHED
                % How many indices matrix has to be built ? (only one if every FUNC=='exp')
                    ind = false(size(indHbH_d)) ;
                    for i = 1:numel(ind) 
                        ind(i) = ~isempty(indHbH_d{i}) ;
                    end
                    ind = ind*[1 0 ; 0 2] ;
                    ind = num2cell(ind,2) ;
                    ind = combvec(ind{:}) ;
                    ind = ind(:,~any(ind==0,1)) ;
                    n_indHbH = size(ind,2) ;
                % Indices matrices
                    indHbH = cell(n_indHbH,1) ;
                    for i = 1:n_indHbH
                        indHbH{i} = indHbH_d{1,ind(1,i)} ;
                        D = 1:length(DIMS_K) ;
                        for d = 2:length(DIMS_K)
                            indHbH{i} = indHbH{i} + (indHbH_d{d,ind(d,i)}-1)*prod(Lk(D<d)) ;
                        end
                    end
            end
        % HBH BUILD
        % Stack The HbH matrices by points of isosurface
            H = zeros(prod(Kk),prod(Mk)*length(indP)) ;
            for p = 1:length(indP)
                sig = Signal(indP(p),:) ;
                SIG = 0 ;
                for i = 1:n_indHbH
                    SIG = SIG + sig(indHbH{i}) ;
                end
                H(:,prod(Mk)*(p-1)+(1:prod(Mk))) = SIG/n_indHbH ;
            end
    end
        
        
% MATRICES F (SHIFT INVARIANCE EVAL.)
    function [F,Wup,Wdwn] = computeF(r,s)
        % R first eigenvectors
            Wp = W(:,1:R0(r)) ;
        % W_up and W_down construction
            shift = SHIFTS(s,:) ;
            D = 1:length(shift) ;
            shift_dims = D(shift~=0) ;
            % Indices Up
                indUp = true(size(indHbH_d{1,1},1),1) ; 
                if any(isCOS.*shift) ; indTrUp = true(size(indHbH_d{1,1},1),1) ; end
                for d = shift_dims
                    switch lower(FUNC(d,:))
                        case 'exp'
                            if shift(d)>0 ; indUp = indUp & (indHbH_d{d,1}(:,1)<=(Kk(d)-shift(d))*DECIM_K(d)) ; end
                            if shift(d)<0 ; indUp = indUp & (indHbH_d{d,1}(:,1)>-shift(d)*DECIM_K(d)) ; end
                        case 'cos'
                            indTrUp = indTrUp & (indHbH_d{d,1}(:,1)<=(Kk(d)-shift(d))*DECIM_K(d)) ;
                            indTrUp = indTrUp & (indHbH_d{d,1}(:,1)>shift(d)*DECIM_K(d)) ;
                    end
                end
            % Partial Matrix Up
                if any(isCOS.*shift) 
                    Wup = Wp(indUp & indTrUp,:) ;
                else
                    Wup = Wp(indUp,:) ;
                end
            % Indices Up
                indDwn = true(size(indHbH_d{1,1},1),1) ; 
                if any(isCOS.*shift) 
                    indTrDwn1 = true(size(indHbH_d{1,1},1),1) ; 
                    indTrDwn2 = true(size(indHbH_d{1,1},1),1) ; 
                end
                for d = shift_dims
                    switch lower(FUNC(d,:))
                        case 'exp'
                            if shift(d)>0 ; indDwn = indDwn & (indHbH_d{d,1}(:,1)>shift(d)*DECIM_K(d)) ; end
                            if shift(d)<0 ; indDwn = indDwn & (indHbH_d{d,1}(:,1)<=(Kk(d)+shift(d))*DECIM_K(d)) ; end
                        case 'cos'
                            indTrDwn1 = indTrDwn1 & (indHbH_d{d,1}(:,1)<=(Kk(d)-2*abs(shift(d)))*DECIM_K(d)) ;
                            indTrDwn2 = indTrDwn2 & (indHbH_d{d,1}(:,1)>2*abs(shift(d))*DECIM_K(d)) ;
                    end
                end
            % Partial Matrix Down
                if any(isCOS.*shift) 
                    Wdwn = (Wp(indDwn & indTrDwn1,:) + Wp(indDwn & indTrDwn2,:))/2 ;
                else
                    Wdwn = Wp(indDwn,:) ;
                end
        % Shift invariance evaluation
            switch FIT
                case 'LS'
                    F = Wup\Wdwn ;
                case 'TLS'
                    [E,~] = svd([Wup' ; Wdwn']*[Wup Wdwn],0) ;
                    F = -E(1:R0(r),R0(r)+(1:R0(r)))/E(R0(r)+(1:R0(r)),R0(r)+(1:R0(r))) ;
            end
    end


% ESTER CRITERION
    function err = ester(r,s)
        [F,Wup,Wdwn] = computeF(r,s) ;
        err = 1/norm(Wup*F-Wdwn) ;
    end


% WAVEVECTORS K EXTRACTION
    function extractPoles(r)
        % Evaluation of all shift invariances
            F = zeros(R0(r),R0(r),size(SHIFTS,1)) ;
            for s = 1:size(SHIFTS,1)
                F(:,:,s) = computeF(r,s) ;
            end
        % Evaluation of PHI
            if size(SHIFTS,1)==1 % One shift, simple diagonalization
                [~,PHI] = eig(F) ;
            else % More than one shift, joint diagonalization
                % Jacobi angles (fails when there is multiplicity)
                    if 0
                        Ft = reshape(F,size(F,1),[]) ;
                        [~,PHI] = joint_diag(Ft,1e-8) ;
                    end
                % Linear combination of F's
                    if 1
                        Beta = 1 + 1.2.^(0:size(SHIFTS,1)-1) ;
                        Beta = Beta/sum(Beta) ;
                        Ft = sum(F.*repmat(reshape(Beta,[1 1 length(Beta)]),[size(F,1) size(F,2) 1]),3) ;
                        [T,~] = eig(Ft) ;
                        PHI = zeros(size(Ft).*[1 size(SHIFTS,1)]) ;
                        for s = 1:size(SHIFTS,1)
                            PHI(:,(s-1)*size(F,2)+(1:size(F,2))) = T\F(:,:,s)*T ;
                        end
                    end
            end
        % Signal Poles
            indDiag = repmat(eye(R0(r)),[1 size(SHIFTS,1)])==1 ;
            Z = reshape(PHI(indDiag),[R0(r) size(SHIFTS,1)]).' ;
        % Wavevectors in the SHIFT basis
            K = log(Z)/1i ;
        % Wavevectors in the cartesian basis
            K = (SHIFTS*diag(DECIM_K))\(K) ;
        % If COSINUS searched...
            K(isCOS,:) = acos(exp(1i*K(isCOS,:))) ;
    end


% AMPLITUDES U RETRIEVAL
    function computeU % Here, decimation factors are not taken into account
        % Vandermonde Matrix
            % Indices
                indV = zeros(prod(Lk),length(Lk)) ;
                D = 1:length(DIMS_K) ;
                for d = 1:length(D)
                    indV(:,d) = kron(repmat((1:Lk(d))',[prod(Lk(D>d)) 1]),ones([prod(Lk(D<d)) 1])) - 1 ;
                end
            % Poles
                Z = exp(1i.*K) ;
                Z = repmat(reshape(Z,[1 size(Z)]),[prod(Lk) 1]) ;
            % Conditionning (using a shift)
                indInv = abs(Z)>1 ;
                shift = 0*indInv.*repmat(1-Lk,[prod(Lk) 1 size(Z,3)]) ;
                indV = repmat(indV,[1 1 size(Z,3)]) + shift ;
            % Vandermonde
                V = reshape(prod(Z.^indV,2),[prod(Lk) size(Z,3)]) ;
        % Signal reshaping
            S = Signal.' ;
        % Amplitude estimation
            A = V\S ;
        % Shift conpensation
            shi = diag(V(1,:)) ;
            A = shi*A ;
        % Reshaping
            U = reshape(A.',[Lp , size(K,2)]) ;
    end
 


    

% ===================================================================================================================    
% FUNCTIONS FOR UTILS
% ===================================================================================================================

% PROCESS INPUTS
    function parseInputs
        nargin = length(varargin) ;
        % Is DIMS_K the first argument ?
            if mod(nargin,2)==1 
                if isnumeric(varargin{1})
                    DIMS_K = varargin{1} ;
                else
                    errorInput('Wrong second argument : should be DIMS_K or a string') ;
                end
                if nargin>1
                    varargin = varargin(2:end) ;
                end
            end
        % Treat following arguments
            paramSet = false(100) ;
            if nargin>1
                for i = 1:2:length(varargin)-1
                    Name = varargin{i} ;
                    Value = varargin{i+1} ;
                    switch upper(Name)
                        case 'DIMS_K'
                            DIMS_K = Value ;
                            paramSet(1) = true ;
                        case 'FUNC'
                            FUNC = Value ;
                            paramSet(2) = true ;
                        case 'R0'
                            R0 = Value ;
                            paramSet(3) = true ;
                        case 'ESTER_THRS'
                            ESTER_THRS = Value ;
                            paramSet(4) = true ;
                        case 'FIT'
                            FIT = upper(Value) ;
                            paramSet(5) = true ;
                        case 'DECIM'
                            DECIM = Value ;
                            paramSet(6) = true ;
                        case 'SHIFTS'
                            SHIFTS = Value ;
                            paramSet(7) = true ;
                        case 'CHOICE'
                            CHOICE = Value ;
                            paramSet(8) = true ;
                        case 'STABILDIAG'
                            STABILDIAG = Value ;
                            paramSet(9) = true ;
                        case 'DEBUG'
                            DEBUG = Value ;
                            paramSet(10) = true ;
                        otherwise
                            errorInput(['Wrong argument name in n°',num2str(i),'.'])
                    end
                end
            end
        % DEFAULT VALUES
            if ~paramSet(1) ; DIMS_K = ndims(Signal) ; end
            if ~paramSet(2) ; FUNC = repmat('exp',[length(DIMS_K) 1]) ; end
            if ~paramSet(3) ; R0 = 1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)) ; end
            if ~paramSet(4) ; ESTER_THRS = 1 ; end
            if ~paramSet(5) ; FIT = 'TLS' ; end
            if ~paramSet(6) ; DECIM = ones(1,ndims(Signal)) ; end
            if ~paramSet(7) ; SHIFTS = eye(length(DIMS_K)) ; end
            if ~paramSet(8) ; CHOICE = 'auto' ; end
            if ~paramSet(9) ; STABILDIAG = false ; end
            if ~paramSet(10) ; DEBUG = false ; end
    end


% PROMPT AN ERROR ON WRONG INPUT ARGUMENTS
    function errorInput(info) 
        msg = ['Input arguments incorrect.',char(10)] ;
        msg = [msg,info] ;
        error([msg,char(10)])
    end


    

% ===================================================================================================================    
% GRAPHICAL FUNCTIONS
% ===================================================================================================================

% PLOT THE STABILIZATION DIAGRAM
    function stabilizationDiagram()
        if length(DIMS_K)>1 ; warning('Stabilization Diagram is available for 1D-ESPRIT only.') ; return ; end
        Kstab = zeros(length(R0),max(R0))*NaN ;
        for r = 1:length(R0)
            extractPoles(r) ; 
            Kstab(r,1:R0(r)) = K ;
        end
        K = Kstab(R,1:R) ;
        Kstab = sort(Kstab,2) ;
        % Figure
            stab.fig = figure('NumberTitle','off','Name','ESPRIT : STABILIZATION DIAGRAM. Click to select an other signal order. Right Click to quit.') ; %figure ;
            % Axes for the signal order selection Criterions
                stab.axCrit = axes('outerposition',[0 0 .2 1]) ;
                    stab.plOrderLine = plot(stab.axCrit,R0(R)*[1 1],[min(log10(ESTER(:))) max(log10(ESTER(:)))],'-.k','linewidth',1) ;
                    stab.plRLine = plot(stab.axCrit,R0(R)*[1 1],[min(log10(ESTER(:))) max(log10(ESTER(:)))],'-.r','linewidth',1) ;
                    plot(R0,log10(ESTER),'.-','markersize',20,'linewidth',1) ;
                    set(gca,'view',[-90 90]) ;
                    box on
                    grid on
                    axis tight
                    % Disable rotate3d
                        hBehavior = hggetbehavior(stab.axCrit, 'Rotate3d');
                        hBehavior.Enable = false ;
            % Axes for the Poles
                stab.axPoles = axes('outerposition',[.2 0 .8 1]) ;
                    plot3(abs(real(Kstab)),repmat(R0(:),[1 max(R0)]),abs(imag(Kstab)./real(Kstab)),'.','markersize',13,'linewidth',1.5)
                    stab.axPoles.ZScale = 'log' ;
                    stab.axPoles.SortMethod = 'childorder' ;
                stab.axMean = axes('position',stab.axPoles.Position) ;
                    meanFFTSignal = Signal ;
                    meanFFTSignal = fft(meanFFTSignal,[],2) ;
                    meanFFTSignal = abs(meanFFTSignal) ;
                    meanFFTSignal = mean(meanFFTSignal,1) ;
                    plot((0:prod(Lk)-1)/prod(Lk)*2*pi,log10(meanFFTSignal),'k') ;
                    axis off
                    box on
            % Link axes
                set([stab.axMean stab.axPoles],'xscale','log')
                set([stab.axMean stab.axPoles],'xlim',[1/prod(Lk)*2*pi/2 pi]) ;
                set([stab.axPoles],'ylim',[R0(1) R0(end)]) ;
                global hlink , hlink = linkprop([stab.axMean stab.axPoles],'position') ;
                linkaxes([stab.axMean stab.axPoles],'x') ;
                uistack(stab.axMean,'bottom') ;
            % Poles Chosen at the end
                stab.plRPoles = plot3(stab.axPoles,abs(real(Kstab(R,:))),R0(R)*ones(1,max(R0)),abs(imag(Kstab(R,:))./real(Kstab(R,:))),'or','linewidth',1.5) ;
            % Poles at the mouse position
                stab.plOrderPoles = plot3(stab.axPoles,abs(real(Kstab(R,:))),R0(R)*ones(1,max(R0)),abs(imag(Kstab(R,:))./real(Kstab(R,:))),'.k','markersize',20) ;
            % Figure Callbacks setting
                stab.fig.WindowButtonMotionFcn = @(src,evt)changeStabDiagOrder(Kstab,stab,'move') ;
                stab.fig.WindowButtonDownFcn = @(src,evt)changeStabDiagOrder(Kstab,stab,'click') ;
            % Buttons for the view
                stab.btnSwitch = uicontrol(stab.fig,'style','pushbutton') ;
                stab.btnSwitch.String = 'Frequency' ;
                stab.btnSwitch.TooltipString = 'Switch Representation Mode' ;
                stab.btnSwitch.Units = 'normalized' ;
                margin = 0.003 ;
                btnWidth = 0.08 ;
                btnHeight = 0.03 ;
                stab.btnSwitch.Position = [1-btnWidth-margin 1-btnHeight-margin btnWidth btnHeight] ;
                stab.btnSwitch.Callback = @(src,evt)btnSwitchCallback(stab) ;
            % Wait for the figure to be closed
                uiwait(stab.fig) ;
    end

% CHANGE THE STABILIZATION DIAGRAM ORDER WITH MOUSE POSITION
    function changeStabDiagOrder(Kstab,stab,event)
        % Get the order given by the mouse position
            mouseOrder = stab.axCrit.CurrentPoint(1,1) ;
            selectOrder = min(max(R0(1),round(mouseOrder)),R0(end)) ;
            [~,order] = min(abs(R0-selectOrder)) ;
            order = order(1) ;
        % Retrieve identified poles
            stab.plOrderLine.XData = R0(order)*[1 1] ;
            stab.plOrderPoles.XData = abs(real(Kstab(order,:))) ;
            stab.plOrderPoles.YData = R0(order)*ones(1,max(R0(:))) ;
            stab.plOrderPoles.ZData = abs(imag(Kstab(order,:))./real(Kstab(order,:))) ;
        % If Clicked, change the final output poles
            switch event
                case 'move'
                case 'click'
                    R = order ;
                    K = Kstab(R,1:R) ;
                    stab.plRLine.XData = R0(order)*[1 1] ;
                    stab.plRPoles.XData = abs(real(Kstab(order,:))) ;
                    stab.plRPoles.YData = R0(order)*ones(1,max(R0(:))) ;
                    stab.plRPoles.ZData = abs(imag(Kstab(order,:))./real(Kstab(order,:))) ;
            end
        % Close the figure if double clicked
            switch stab.fig.SelectionType
                case 'open'
                case 'alt'
                    close(stab.fig) ;
            end
    end

% CHANGE STABIL. DIAGRAM ROTATION
    function btnSwitchCallback(stab)
        switch stab.btnSwitch.String
            case 'Frequency' % switch to Damping
                stab.axPoles.View = [0 0] ;
                stab.btnSwitch.String = 'Damping' ;
            case 'Damping' % switch to Frequency
                stab.axPoles.View = [0 90] ;
                stab.btnSwitch.String = 'Frequency' ;
        end
    end

end