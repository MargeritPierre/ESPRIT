function OUT = ESPRIT(Signal,varargin)
%
% ESPRIT Extract wavevectors from a signal
%
%   OUT = ESPRIT(Signal,varargin)
%
%   varargin parameters (and default values) : [] set the parameter to default
%       CHOICE : not implemented 
%           string : 'auto' or 'manual' ('auto')
%       SOLVER : signal subspace estimation
%           string : 'eigs' or 'eig' ('eig') 
%       CRITERION : signal order choice criterion
%           string : 'ESTER', 'SAMOS', 'MDL' or 'ALL' ('MDL')
%       CRIT_THRS : signal order criterion threshold 
%           1x1 double (1)
%       COMPUTE_U : compute amplitudes 
%           boolean (false)
%       SIGNAL_MODEL : compute the signal model 
%           boolean (false)
%       COMPUTE_dK : compute uncertainties on wavevctors
%           boolean (false)
%       COMPUTE_dU : compute uncertainties on amplitudes 
%           boolean (false)
%       DIMS_K : dimensions of the wavevectors 
%           1xD integer (ndims(Signal))
%       FUNC : type of function searched 
%           Dx3 char : 'exp' or 'cos' (repmat('exp',[length(DIMS_K) 1]))
%       R0 : signal order candidates 
%           1xnR integer : (1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)))
%       M/L : Spatial Smoothing ratio (gives the HbH matrix shape)
%           1xD double : 0<M_d/L_d<1 (2/3)
%               - 0 : no spatial smoothing (M = ones(1,D))
%               - scalar : same M/L ratio for all dimensions
%               - 1xD double : different for each dimension 
%       W0 : signal subspace of a previous iteration 
%           Prod(K)xMax(R0) double : ([])
%               - [] for a random subspace initialization
%       FIT : regression method for the spectral matrix estimation
%           string : 'LS' or 'TLS' (TLS)
%       DECIM : decimation factors 
%           1xD integer : (ones(1,ndims(Signal)))
%           - along DIMS_K : /!\ NYQUIST
%           - along DIMS_P : reduce the size of Css.
%           Do not affect the estimation of U.
%       SHIFTS : for multiresolution 
%           QxD integer : (eye(length(DIMS_K)))
%       DEBUG : prompt procedure state 
%           boolean (false)
%       STABILDIAG : plot the stabilization diagram 
%           boolean (false)
%       MAC : correlation coefficient to link the poles by MAC values 
%           1x1 double : (0)
%           - MAC == 0 : no MAC value computed
%           - 0 < MAC < 1 : poles with corr > MAC are linked
%       Kstab : only show the stabil. diag. for previous results
%           length(R0)xmax(R0) complex double with NaNs ([])
%       Ustab : only show the stabil. diag. for previous results
%           prod(Lp)xlength(R0)xmax(R0) complex double with NaNs ([])
%       STOP : partial execution of the algorithm
%           string ('None') : 
%           - 'Hankel'
%           - 'Covariance'
%           - 'Subspace'
%           - 'Order'
%           - 'Wavevectors'
%           - 'Amplitudes'
%           - 'None'
%
%
%   OUT structure : 
%    % Sizes
%       Lp : dimensions of the signal over the isophase surface (size : [1 length(DIMS_P)])
%       Lk : dimensions of the signal over the wavevectors (size : [1 length(DIMS_K)])
%       Kk : number of rows of the signal matrix for each signal dimension (size : [1 length(DIMS_K)])
%       Mk : number of columns of the signal matrix for each signal dimension (size : [1 length(DIMS_K)])
%    % HbH Matrix
%       indHbH : cell of indices for the HbH matrix
%       indP : indices of the isophase surface points used to compute the covariance matrix Css
%    % Subspace
%       Css : Autocovariance matrix (size : [prod(Kk) prod(Mk)])
%       W : Truncated signal subspace, max(R0) first E.Vec of CSS (size : [prod(Kk) max(R0)])
%       lambda : max(R0) first E.Val of Css (size : [max(R0) 1])
%    % Signal Order Criterion
%       CRIT : final signal order criterion (size : [1 length(R0)])
%       ESTER  : ester criterion for each shift (size : [size(SHIFTS,1) length(R0)])
%    % Results
%       K : extracted complex wavevectors (size : [length(DIMS_K) R])
%       dK : uncertainties (size : [length(DIMS_K) R])
%       U : complex amplitudes (modes) (size : [[Lp] R])
%
%
% This ESPRIT algorithm implements :
%   - Standard ('exp') and Paired ('cos') ESPRIT (FUNC)
%   - ESTER criterion (ESTER_THRS)
%   - Multidimensionnal (DIMS_K)
%   - LS or TLS linear regression (FIT)
%   - Multiple Invariance or Multiple Resolution (SHIFTS)
%   - Decimate ESPRIT (DECIM)
%   - Stabilization Diagram (STABILDIAG)
%   - Choice between EIG, EIGS computing W for speed
%   - Allow for the signal subspace tracking for speed (trough W and W0)
%    
% To go further on developments :
%   (1) use of high-order cummulents
%   (2) MR-ESPRIT : solving ambiguity when Nyquist is violated
%   (3) GUI for manual choice
%   (4) Cramer-Rao bound, uncertainty estimation
%   (5) replace kron() by repelem() in indHbH_d build for speed (ver>=Matlab2015)


% INITIALIZATION
        
   % Output initialization
        OUT = [] ;
        OUT.Title = 'ESPRIT_Results' ;
        OUT.Signal = Signal ;

    % Input Initialization
        parseInputs ; 
        varargin ;
        paramSet ;
        CHOICE ;
        SOLVER ;
        CRITERION ;
        CRIT_THRS ; 
        COMPUTE_U ;
        COMPUTE_dK ;
        DIMS_K ; 
        FUNC ; 
        R0 ; 
        M_L ;
        W0 ;
        FIT ; 
        DECIM ; 
        SHIFTS ; 
        DEBUG ;
        STABILDIAG ;
        MAC ;
        Kstab ;
        Ustab ;

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
        NDIMS = ndims(Signal) ; DIMS = 1:NDIMS ; % Dimension of the signal
        DIMS_P = DIMS(~ismember(DIMS,DIMS_K)) ; % Dimensions of the isosurface (avoids setdiff() which is slow)
        Lk = SIZE(DIMS_K) ; % Number of points in the dimensions of the exponentials
        Lp = SIZE(DIMS_P) ; % Number of points in the dimensions of the isophase surface
        % When only one point over the isophase surface
            if isempty(DIMS_P) ; Lp = 1 ; end
        % Decimation
            DECIM_K = DECIM(DIMS_K) ;
            DECIM_P = DECIM(DIMS_P) ;
        % Directions of exponentials & cosinuses
            isCOS = ismember(num2cell(FUNC,2)','cos') ;
            isEXP = ~isCOS ;
            
    % Reshaping the signal
        Signal = permute(Signal,[DIMS_P DIMS_K]) ; % Sort dimensions
        Signal = reshape(Signal,[prod(Lp) prod(Lk)]) ;
        [~,permDims] = sort([DIMS_P DIMS_K]) ; % Used to rechape the signal model
    
    % Mean Value
        %Signal = Signal-repmat(mean(Signal,2),[1 size(Signal,2)]) ;
        
        
    % Output Initialization
        K = [] ; % Wavevectors Matrix
        dK = [] ; % Wavevectors Uncertainty
        V = [] ; % Steering Matrix
        U = [] ; % Amplitudes
        dU = [] ; % Amplitudes Uncertainty
        if any(isCOS) ; PHI = [] ; end % Phases of the cosinuses
        SignalModel = [] ; % Reconstructed Signal Model
        dS = [] ; % Errors
        Hss = [] ;
        
        
% ONLY PLOT SOME PREVIOUS RESULTS ON THE STABILIZATION DIAGRAM ?
    if ~isempty(Kstab)
        onlyStabDiag = true ;
        stabilizationDiagram() ;
        return ;
    else
        onlyStabDiag = false ;
    end
        
        
% HANKEL MATRIX BUILDING
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Hankel Matrix Indices : ') ; lastTime = tic ; end
    
    % Determine the hankel matrix sizes
        Kk = [] ;
        Mk = [] ;
        setHankelSize() ;
        
    % Buil sub-Hankel matrices indices
        n_indHbH = any(isCOS)+1 ;
        subIndH = {} ;
        subHankelInd() ;
        
    % Hankel-block-Hankel matrix indices
        
        % Indices for each dimension (used in shift invariances and Vandermonde matrix)
            indHbH_d =  {} ;
            indicesHbH_d() ;
            
        % Complete indices
            indHbH = [] ;
            indP = [] ;
            indicesHbH() ;
            
        
        
        
% AUTOCOVARIANCE MATRIX
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Computing of Css : ') ; lastTime = tic ; end
        % Computation speedup
            if prod(Kk)<prod(Mk)*length(indP)
                transpCss = false ; % Css = H*H' , [W,lambda] = eig(Css)
                Css = buildCss() ;
            else
                transpCss = true ; % Css = H'*H , [Wr,lambda] = eig(Css), W = H*W\diag(sqrt(lambda)) ,
                if(DEBUG) ; display('        transpCss: build Hss') ; end
                Hss = buildHss(Signal) ;
                if(DEBUG) ; display('         ... and Css = Hss''*Hss') ; end
                Css = Hss'*Hss ;
            end
        
        
        
% SIGNAL SUBSPACE ESTIMATION        
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Eigendecomposition of Css : ') ; lastTime = tic ; end
    % Adjust signal order candidates
        R0 = R0(R0<=size(Css,1)-1-max(sum(abs(SHIFTS),2))) ;
    % ESTIMATION
        if isempty(W0) % re-compute the signal subspace
            switch SOLVER
                case 'eig'
                    [W,lambda] = eig(Css,'vector') ;
                case 'eigs'
                    [W,lambda] = eigs(Css,max(R0),'lm') ;
                    lambda = diag(lambda) ;
            end
        else % Approximate with one QR iteration (Badeau-style)
            Cxy = Css*W0 ;
            [W,~] = qr(Cxy,0) ;
            lambda = diag(W'*Css*W) ;
        end
    % SORT EIGENVALUES IN DESCENDING ORDER
        [~,Ind] = sort(lambda,'descend') ;
        lambda = lambda(Ind) ;
        W = W(:,Ind) ;
    % GET THE LEFT EIGENVECTORS
        if transpCss
            if(DEBUG) ; display('        transpCss: compute left eig.vec') ; end
            W = Hss*W ;
            W = bsxfun(@(w,l)w./l, W, sqrt(lambda(:)).') ;
        end
    % KEEP the needed E.V only
        lambda = lambda(1:max(R0)) ;
        W = W(:,1:max(R0)) ;
        
        
% SIGNAL ORGER CRITERION
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Signal Order Criterion : ') ; lastTime = tic ; end
    if length(R0)==1 % THE SIGNAL ORDER HAS BEEN GIVEN
        R = 1 ;
    else % SIGNAL ORDER CHOICE
        % INITIALIZATION
            ESTER = NaN*ones(size(SHIFTS,1),length(R0)) ;
            SAMOS = NaN*ones(size(SHIFTS,1),length(R0)) ;
        % MDL
            MDL = (lambda(1:end-1)-lambda(2:end))./lambda(1:end-1) ;
            MDL = [MDL ; 0] ;
            MDL = MDL(R0) ;
            MDL = MDL(:)' ;
            CRIT = MDL ;
        % OTHER CRITERIONS
            % ESTER
                if strcmp(CRITERION,'ALL') || strcmp(CRITERION,'ESTER')
                    if(DEBUG) ; display('        ester') ; end
                    t = tic ;
                    for r = 1:length(R0)
                        for s = 1:size(SHIFTS,1)
                            ESTER(s,r) = ester(r,s) ;
                        end
                        if length(R0)>1 && r==1 && DEBUG ; disp(['            estim: ',num2str(toc(t)*length(R0)),' secs']) ; end
                    end
                end
            % SAMOS
                if strcmp(CRITERION,'ALL') || strcmp(CRITERION,'SAMOS')
                    if(DEBUG) ; display('        samos') ; end
                    t = tic ;
                    for r = 1:length(R0)
                        for s = 1:size(SHIFTS,1)
                            SAMOS(s,r) = samos(r,s) ;
                        end
                        if length(R0)>1 && r==1 && DEBUG ; disp(['            estim: ',num2str(toc(t)*length(R0)),' secs']) ; end
                    end
                end
        % CHOOSE THE SIGNAL ORDER
            CRIT = prod(CRIT,1) ;
            R = find(CRIT>=max(CRIT)*CRIT_THRS,1,'last') ;
            OUT.CRIT = CRIT ;
    end
    
        
% STABILIZATION DIAGRAM
    if STABILDIAG && length(R0)>1
        Kstab = [] ;
        Ustab = [] ;
        MAC ;
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Stabilization Diagram : ') ; lastTime = tic ; end
        stabilizationDiagram() ;
    end
        
    
% EXTRACT THE WAVEVECTORS
    if ~STABILDIAG || COMPUTE_dK
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Pole Estimation : ') ; lastTime = tic ; end
        [T,~] = extractPoles(R) ;
    end
        
    
% ESTIMATE AMPLITUDES
    if COMPUTE_U
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Amplitude Estimation : ') ; lastTime = tic ; end
        computeU() ;
        A ;
    end

% BUILD THE SIGNAl MODEL
    if SIGNAL_MODEL
        computeSignalModel() ;
    end
    
    
% ESTIMATE UNCERTAINTIES
    if COMPUTE_dK || COMPUTE_dU 
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Uncertainty Estimation : ') ; lastTime = tic ; end
        computeUncertainties ;
    end
        
    
% OUTPUT
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Output : ') ; lastTime = tic ; end
    % Input arguments
        OUT.varargin = varargin ;
        OUT.R0 = R0 ;
    % Results
        OUT.K = K ;
        OUT.V = V ;
        OUT.U = U ;
        if any(isCOS) ; OUT.PHI = PHI ; end
        if ~isempty(SignalModel)
            if isempty(DIMS_P)
                OUT.SignalModel = reshape(SignalModel,Lk) ;
                OUT.dS = reshape(dS,Lk) ;
            else
                OUT.SignalModel = permute(reshape(SignalModel,[Lp Lk]),permDims) ;
                OUT.dS = permute(reshape(dS,[Lp Lk]),permDims) ;
            end
            OUT.RelativeError = norm(dS(:))/norm(Signal(:)) ;
        end
    % Uncertainties
        OUT.dK = dK ;
        OUT.dU = dU ;
    % Sizes
        OUT.Lp = Lp ;
        OUT.Lk = Lk ;
        OUT.Kk = Kk ;
        OUT.Mk = Mk ;
        OUT.isCOS = isCOS ;
        OUT.isEXP = isEXP ;
    % Decimation
        OUT.DECIM_K = DECIM_K ;
        OUT.DECIM_P = DECIM_P ;
    % HbH indices
        OUT.subIndH = subIndH ;
        OUT.indHbH = indHbH ;
        OUT.indP = indP ;
    % Subspace
        OUT.transpCss = transpCss ;
        OUT.Css = Css ;
        OUT.lambda = lambda ;
        OUT.W = W ;
    % Signal Order
        if length(R0)>1
            OUT.CRIT = CRIT ;
            OUT.MDL = MDL ;
            OUT.ESTER = ESTER ;
            OUT.SAMOS = SAMOS ;
        end
    % Stabilization Diagram 
        if STABILDIAG && length(R0)>1
            OUT.Kstab = Kstab ;
            OUT.Ustab = Ustab ;
        end
    

% END OF THE PROCEDURE
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('---------------------------------') ; end

    
    
    
    
    
% ===================================================================================================================    
% FUNCTIONS FOR ESPRIT
% ===================================================================================================================

% HANKEL MATRIX SIZES
    function setHankelSize()
        % Spatial smoothing ratio M/L
            if numel(M_L)==1 ; M_L = M_L*ones(1,length(DIMS_K)) ; end
            M_L = max(M_L,1./Lk) ; % avoids M=0
            M_L = min(M_L,1-1./Lk) ; % avoids K=0
            %M_L = max(M_L,1-M_L) ; % uses the symmetry to M=L/2 to minimize K
        % Sub-grids size
            Kk = floor((1-M_L).*Lk./DECIM_K + 1) ; zeros([1 length(DIMS_K)]) ; % max size of (signal+noise) subspace
        % When Decimation, force the spatial smoothing to use all data
            decimAndNoSS = Lk-(Kk-1).*DECIM_K<DECIM_K ;
            Kk(decimAndNoSS) = Kk(decimAndNoSS)-1 ;
        Mk = Lk - DECIM_K.*(Kk-1) ;
        Kk(isCOS) = floor(Kk(isCOS)/2) ;
    end

% SUB-HANKEL MATRICES INDICES
    function subHankelInd()
        subIndH = cell(length(DIMS_K),n_indHbH) ;
        for d = 1:length(DIMS_K)
            switch lower(FUNC(d,:))
                case 'exp'
                    subIndH{d,1} = (1:DECIM_K(d):DECIM_K(d)*(Kk(d)-1)+1)'*ones(1,Mk(d)) + ones(Kk(d),1)*(0:Mk(d)-1) ;
                case 'cos'
                    subIndH{d,1} = ((Kk(d)-1)*DECIM_K(d)+1:-DECIM_K(d):1)'*ones(1,Mk(d)) + ones(Kk(d),1)*(0:Mk(d)-1) ;
                    subIndH{d,2} = (Kk(d)*DECIM_K(d):DECIM_K(d):DECIM_K(d)*(2*Kk(d)-1))'*ones(1,Mk(d)) + ones(Kk(d),1)*(1:Mk(d)) ;
            end
        end
    end

% HbH INDICES FOR EACH DIMENSION
    function indicesHbH_d()
        indHbH_d = cell(length(DIMS_K),n_indHbH) ;
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
    end

% HANKEL-BLOCK-HANKEL MATRIX INDICES
    function indicesHbH
        % Decimation in the isosurface points
            indP = 1:size(Signal,1) ;
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
                    indHbH = {indHbH_d{1,:}} ;
                    D = 1:length(DIMS_K) ;
                    for d = 2:length(DIMS_K)
                        indHbH{1} = indHbH{1} + (indHbH_d{d,1}-1)*prod(Lk(D<d)) ;
                        if n_indHbH==2
                            if isCOS(d)
                                indHbH{2} = indHbH{2} + (indHbH_d{d,2}-1)*prod(Lk(D<d)) ;
                            else
                                indHbH{2} = indHbH{2} + (indHbH_d{d,1}-1)*prod(Lk(D<d)) ;
                            end
                        end
                    end
            end
    end

% BUILD THE CSS MATRIX
    function Css = buildCss()
        % Build the covariance matrix
            if any(Kk~=Lk) % Mk~=1, spatial smoothing
                if 0 %length(Kk)==1 && sum(isCOS)==0 % 1D exp case, convolution
                    % Signal vector
                        sig = Signal(indP,indHbH{1}(:)) ;
                        for i = 2:n_indHbH
                            sig = sig + Signal(indP,indHbH{i}(:)) ;
                        end
                        sig = sig/n_indHbH ;
                        sig = sig.' ; % now size(sig)=[Kk*Mk length(indP)]
                    % Signal Matrix
                        SIG = reshape(sig,[Kk Mk length(indP)]) ;
                else % N-D Case
                    % Sum the elementary Css matrices by points of isosurface
                        Css = zeros(prod(Kk),prod(Kk)) ;
                        t = tic ; 
                        for p = 1:length(indP)
                            sig = Signal(indP(p),:) ;
                            SIG = sig(indHbH{1}) ;
                            for i = 2:n_indHbH
                                SIG = SIG + sig(indHbH{i}) ;
                            end
                            SIG = SIG/n_indHbH ;
                            Css = Css + SIG*SIG' ;
                            if length(indP)>1 && p==1 && DEBUG ; disp(['       estim: ',num2str(toc(t)*length(indP)),' secs']) ; end
                        end
                end
            else % Mk==1, no spatial smoothing (multiple snapshots case)
                SIG = Signal.' ;
                Css = SIG*SIG' ;
            end
        % Normalize
            %Css = Css/(prod(Mk)*length(indP)) ;
    end


% BUILD THE HSS MATRIX
    function Hss = buildHss(Data)
            if any(Kk~=Lk) % Mk~=1, spatial smoothing
                % Stack The HbH matrices by points of isosurface
                    Hss = zeros(prod(Kk),prod(Mk)*length(indP)) ;
                    for p = 1:length(indP)
                        sig = Data(indP(p),:) ;
                        SIG = sig(indHbH{1}) ;
                        for i = 2:n_indHbH
                            SIG = SIG + sig(indHbH{i}) ;
                        end
                        SIG = SIG/n_indHbH ;
                        Hss(1:prod(Kk),(1:prod(Mk))+(p-1)*prod(Mk)) = SIG ;
                    end
            else % Mk==1, no spatial smoothing (multiple snapshots case)
                Hss = Data(indP,:).' ;
            end
    end


% SELECTION INDICES (AND MATRICES OPTIONALLY)
    function [indUp,indDwn1,indDwn2,Jup,Jdwn] = selectMatrices(shift)
            % Infos
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
            % Indices Down
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
            % Combine Indices
                if any(isCOS.*shift) % Cosinus searched in the shift
                    indUp = indUp & indTrUp ;
                    indDwn1 = indDwn & indTrDwn1 ;
                    indDwn2 = indDwn & indTrDwn2 ;
                else % Only exponentials in the shift
                    indDwn1 = indDwn ;
                    indDwn2 = [] ;
                end
            % Build selection MATRICES
                if nargout>3 % If Jup and Jdwn are needed
                    I = speye(prod(Kk)) ;
                    % Jup
                        Jup = I(indUp,:) ;
                    % Jdwn
                        Jdwn = I(indDwn1,:) ;
                        if any(isCOS.*shift) % Cosinus searched in the shift
                            Jdwn2 = I(indDwn2,:) ;
                            Jdwn = (Jdwn + Jdwn2)/2 ;
                        end
                end
    end
        
        
% MATRICES F (SHIFT INVARIANCE EVAL.)
    function [Wup,Wdwn,F] = computeF(r,s)
        % Get indices
            shift = SHIFTS(s,:) ;
            [indUp,indDwn1,indDwn2] = selectMatrices(shift) ;
        % R first eigenvectors
            Wp = W(:,1:R0(r)) ;
        % W_up and W_down construction
            if any(isCOS.*shift) 
                Wup = Wp(indUp,:) ;
                Wdwn = (Wp(indDwn1,:) + Wp(indDwn2,:))/2 ;
            else
                Wup = Wp(indUp,:) ;
                Wdwn = Wp(indDwn1,:) ;
            end
        % Only the Wup and Wdwn matrices needed ?
            if nargout==2 ; return ; end
        % Shift invariance evaluation
            switch FIT
                case 'LS'
                    F = Wup\Wdwn ;
                case 'TLS'
                    [E,~] = svd([Wup' ; Wdwn']*[Wup Wdwn],0) ;
                    F = -E(1:R0(r),R0(r)+(1:R0(r)))/E(R0(r)+(1:R0(r)),R0(r)+(1:R0(r))) ;
                case 'TLS2'
                    L_2 = floor(Lk./2) ;
                    Wup = Wp(1:2:L_2*2,:) ;
                    Wdwn = Wp(2:2:L_2*2,:) ;
                    [E,~] = svd([Wup' ; Wdwn']*[Wup Wdwn],0) ;
                    F = -E(1:R0(r),R0(r)+(1:R0(r)))/E(R0(r)+(1:R0(r)),R0(r)+(1:R0(r))) ;
            end
    end


% ESTER CRITERION
    function err = ester(r,s)
        [Wup,Wdwn,F] = computeF(r,s) ;
        err = 1/norm(Wup*F-Wdwn) ;
    end


% SAMOS CRITERION
    function err = samos(r,s)
        [Wup,Wdwn] = computeF(r,s) ;
        S = [Wup Wdwn] ;
        g = real(sort(sqrt(eig(S'*S,'vector')))) ;
        %g = flip(svd(S)) ;
        E = sum(g(1:R0(r)))/R0(r) ;
        err = 1/E ;
    end


% WAVEVECTORS K EXTRACTION
    function [T,Beta] = extractPoles(r)
        % Evaluation of all shift invariances
            F = zeros(R0(r),R0(r),size(SHIFTS,1)) ;
            for s = 1:size(SHIFTS,1)
                [~,~,F(:,:,s)] = computeF(r,s) ;
            end
        % Evaluation of PHI
            if size(SHIFTS,1)==1 % One shift, simple diagonalization
                [T,PI] = eig(F) ;
                Beta = 1 ;
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
                        PI = zeros(size(Ft).*[1 size(SHIFTS,1)]) ;
                        for s = 1:size(SHIFTS,1)
                            PI(:,(s-1)*size(F,2)+(1:size(F,2))) = T\F(:,:,s)*T ;
                        end
                    end
            end
        % Signal Poles
            indDiag = repmat(eye(R0(r)),[1 size(SHIFTS,1)])==1 ;
            Z = reshape(PI(indDiag),[R0(r) size(SHIFTS,1)]).' ;
        % Wavevectors in the SHIFT basis
            shiftsCOS = any(logical(repmat(isCOS(:)',[size(SHIFTS,1) 1])).*SHIFTS,2) ;
            K = zeros(size(Z)) ;
            K(~shiftsCOS,:) = log(Z(~shiftsCOS,:))/1i ; % FUNC = 'EXP' ;
            K(shiftsCOS,:) = acos(Z(shiftsCOS,:)) ; % FUNC = 'COS' ;
        % Wavevectors in the cartesian basis
            K = (SHIFTS*diag(DECIM_K))\(K) ;
    end



% VANDERMONDE MATRIX
    function buildVandermonde(L) % Here, decimation factors are no longer taken into account
        % Indices
            indV = zeros(prod(L),length(L)) ;
            D = 1:length(DIMS_K) ;
            for d = 1:length(D)
                indV(:,d) = kron(repmat((1:L(d))',[prod(L(D>d)) 1]),ones([prod(L(D<d)) 1])) - 1 ;
            end
        % Wavevectors
            if ~any(isCOS) % No cosinuses, normal Vandermonde matrix
                Kt = K ;
            else % Cosinuses searched, extended Vnadermonde matrix
                Kt = [K bsxfun(@times,double(~isCOS(:)),K)-bsxfun(@times,double(isCOS(:)),K)] ;
            end
        % Vandermonde
            Kt = permute(Kt,[3 2 1]) ;
            indV = permute(indV,[1 3 2]) ;
            V = exp(1i*sum(bsxfun(@times,indV,Kt),3)) ;
    end


% AMPLITUDES U RETRIEVAL
    function computeU % Here, decimation factors are no longer taken into account
        % Vandermonde Matrix
            if isempty(V) ; buildVandermonde(Lk) ; end
        % Shift for the conditionning
            %v0 = max(abs(V),[],1) ;
            %V = V*diag(1./v0) ;
        % Signal reshaping
            S = Signal.' ;
        % Amplitude estimation
            A = V\S ;
        % Vandermonde Shift compensation
            %A = diag(1./v0)*A ;
        % Reshaping
            if ~any(isCOS) % No cosinuses, no phases to estimate
                U = reshape(A.',[Lp , size(K,2)]) ;
            else
                Bplus = A(1:end/2,:) ;
                Bminus = A(end/2+1:end,:) ;
                U = reshape(2*sqrt(Bplus.*Bminus).',[Lp , size(K,2)]) ;
                PHI = reshape(1i/2*log(Bminus./Bplus).',[Lp , size(K,2)]) ;
            end
    end

% SIGNAL MODEL CONSTRUCTION
    function computeSignalModel()
        if isempty(V) ; computeU() ; end
        SignalModel = (V*A).' ;
        dS = Signal-SignalModel ;
    end


% UNCERTAINTY ESTIMATION dK
    function computeUncertainties
        % OPTIONS (hard-coded for now)
            hugeData = true ; % handles huge Data by dividing the data at points
            estimate =  ... 'std' ... % Standard Deviation
                         'delta' ... % Sensibility or perturbation
                        ;
            lin_method =     'conv' ... % Linearization method for the variance:
                            ... 'kron'...
                            ... 'none'...  
                        ;
            covar_estim =   ... 'uniform' ... % Estimation of the perturbation covariance: 
                             'diagonal' ...
                            ;
            formulation =   ... 'analytic' ... % with Vandermonde matrices etc. Do not work with cosinuses
                             'eigenspace' ... % from the eigenspace
                            ;
        % SIGNAL PERTURBATION AND COVARIANCE
            if isempty(SignalModel) ; computeSignalModel() ; end
            if ~strcmp(lin_method,'none') && strcmp(estimate,'std') % standard deviation estimation: std
                if(DEBUG) ; display('        data covariance') ; end
                switch covar_estim
                    case 'uniform'
                        var_dS = var(dS(:),0); % scalar
                    case 'diagonal'
                        var_dS = abs(dS).^2 ; % matrix !
                end
            end
        % WAVEVECTOR PERTURBATION
            if COMPUTE_dK
                if(DEBUG) ; display('        wavevector uncertainties') ; end
                % Perturbed HbH matrix if needed
                    if strcmp(lin_method,'none') %|| strcmp(estimate,'delta')
                        if(DEBUG) ; display('             dHss') ; end
                        dH = buildHss(dS) ;
                    end
                % Selection matrices if needed
                    if(DEBUG) ; display('             selection matrices') ; end
                    switch lin_method
                        case 'kron'
                            % M Matrix
                                Ip = speye(length(indP)) ;
                                I = speye(prod(Lk)) ;
                                M = sparse(prod(Kk)*prod(Mk)*length(indP),prod(Lk)*length(indP)) ;
                                for m = 1:prod(Mk)
                                    Im = I(indHbH{1}(:,m),:) ;
                                    for i = 2:n_indHbH
                                        Im = Im + I(indHbH{i}(:,m),:) ;
                                    end
                                    Mm = kron(Ip,Im);
                                    M((1:prod(Kk)*length(indP))+(m-1)*(prod(Kk)*length(indP)),:) = Mm ;
                                end
                        otherwise
                            % Jdelta matrix (decimation selection matrix)
                                I = speye(prod(Lk),prod(Lk-Mk+1)) ;
                                Jdelta = I(indHbH{1}(:,1),:) ;
                                for i = 2:n_indHbH
                                    Jdelta = Jdelta/2 + I(indHbH{i}(:,1),:)/2 ;
                                end
                    end
                % INITIALIZATION
                    dK = zeros(size(K)) ;
                    shifts = eye(length(DIMS_K)) ; % /!\ SHIFTS IS DEFAULT HERE !
                % HUGE DATA LOOP
                    if hugeData 
                        indP_bkp = indP(:)' ; % the for-loop will iterate
                    else
                        indP_bkp = indP(:) ; % the for-loop will NOT iterate
                    end
                % GOOOOOO !
                    for indP = indP_bkp
                        % Pre-Computations
                            % Right vector of the bilinear form
                                if(DEBUG) ; display('             right eigenvectors of Css') ; end
                                switch formulation
                                    case 'analytic'
                                        % Partial Vandermonde matrices
                                            if isempty(V) ; buildVandermonde(Lk) ; end
                                            P = V(indHbH{1}(:,1),:) ;
                                            Q = V(indHbH{1}(1,:),:) ;
                                        % Complete right-Vandermonde Matrix
                                            QA = zeros(prod(Mk)*length(indP),R0(R)) ;
                                            for p = 1:length(indP) 
                                                QA((1:prod(Mk))+(p-1)*prod(Mk),:) = Q*diag(A(:,indP(p))) ;
                                            end
                                            x = (QA\eye(size(QA,1)))' ;
                                    case 'eigenspace'
                                        if isempty(Hss) ; Hss = buildHss(SignalModel) ; end
                                        x = conj(bsxfun(@(x,c)x./c,Hss'*W(:,1:R0(R)),lambda(1:R0(R)).')*T) ;
                                end
                        % Loop over the shifts
                            if(DEBUG) ; display('             uncertainties') ; end
                            for s = 1:size(K,1)
                                [~,~,~,Jup,Jdwn] = selectMatrices(shifts(s,:)) ;
                                % pre-compute
                                    switch formulation
                                        case 'analytic'
                                            vn = (Jup*P)\speye(size(Jdwn,1)) ;
                                        case 'eigenspace'
                                            vn = T\((Jup*W(:,1:R0(R)))\speye(size(Jdwn,1))) ;
                                    end
                                % Loop over the R0(R) components
                                    for r = 1:size(K,2)
                                        % The polar component
                                            switch FUNC(s,:)
                                                case 'exp'
                                                    PIrn = exp(1i*K(s,r)*DECIM_K(s)) ;
                                                    Arn = exp(1i*K(s,r)*DECIM_K(s)) ;
                                                case 'cos'
                                                    PIrn = cos(K(s,r)*DECIM_K(s)) ;
                                                    Arn = abs(sin(K(s,r)*DECIM_K(s))) ;
                                            end
                                        % Linearization method
                                            switch lin_method
                                                case 'none' % NO LINEARIZATION, uncertainty only
                                                    vrn = (vn(r,:)*(Jdwn-PIrn*Jup))' ;
                                                    dK(s,r) = dK(s,r) + abs(vrn.'*dH*x(:,r)) ;
                                                case 'conv' % LINEAR / BY CONVOLUTION (USE OF THE HANKEL SHAPE OF Hss)
                                                    vrn = (vn(r,:)*(Jdwn-PIrn*Jup)*Jdelta)' ;
                                                    VRN = full(reshape(vrn,[1 Lk-Mk+1])) ;
                                                    XR = reshape(x(:,r),[length(indP) Mk]) ;
                                                    ZN = ifftn(bsxfun(@times,fftn(VRN,[1 Lk]),fftn(XR,[length(indP) Lk]))) ; % ND convolution
                                                    zn = conj(ZN(:)) ;
                                                case 'kron' % LINEAR / BY VECTORIZATION (USES vec(A*X*B) = kron(B.',A)*vec(X) ) 
                                                    vrn = (vn(r,:)*(Jdwn-PIrn*Jup))' ;
                                                    zn = (kron((x(:,r)),vrn)'*M).' ;
                                            end
                                        % Perturbation estimate
                                            if ~strcmp(lin_method,'none')
                                                switch estimate
                                                    case 'std'
                                                        switch covar_estim
                                                            case 'uniform'
                                                                dK(s,r) = dK(s,r) + var_dS*sum(abs(zn).^2) ;
                                                            case 'diagonal'
                                                                dK(s,r) = dK(s,r) + sum(abs(zn).^2.*reshape(var_dS(indP,:),[],1)) ;
                                                        end
                                                    case 'delta'
                                                        dK(s,r) = dK(s,r) + abs(zn)'*reshape(dS(indP,:),[],1) ;
                                                end
                                            end  ;
                                    end % end of this signal order r
                            end % end of this shift s
                    end % end of the points indP
                % Common terms
                    if strcmp(estimate,'std') ; dK = sqrt(dK) ; end
                    dK = diag(DECIM_K)*abs(dK)/abs(Arn) ;
                % Backup the point indices
                    indP = indP_bkp(:)' ;
            end
        % AMPLITUDES PERTURBATION
            if COMPUTE_dU
                if length(DIMS_K)>1 || any(isCOS(:)) ; dU = NaN*ones(size(U)) ; return ; end
                if(DEBUG) ; display('        amplitudes uncertainty') ; end
                if isempty(V) ; buildVandermonde(Lk) ; end
                % Pre-compute the inverse of V
                    if(DEBUG) ; display('             inverse of V') ; end
                    invV = V\eye(prod(Lk)) ;
                % Uncertainty
                    switch estimate
                        case 'std'
                            dU = zeros(R0(R),prod(Lp));
                            for r = 1:R0(R)
                                zrn = invV(r,:)' ;
                                    switch covar_estim
                                        case 'uniform'
                                            dU(r,:) = sqrt(var_dS*sum(abs(zrn).^2)) ;
                                        case 'diagonal'
                                            dU(r,:) = sqrt(var_dS*abs(zrn(:)).^2).' ;
                                    end
                                dU = dU*2 ; % Yes, I don't know why this factor 2...
                            end
                        case 'delta'
                            dV = bsxfun(@times,(0:Lk-1)',V)*diag(dK(1,:)) ;
                            dU = abs(abs(invV)*(dS.' + abs(dV)*U.'));
                    end
                % Final Processing
                    dU = reshape(dU.',[Lp , size(K,2)]) ;
            end
    end


    

% ===================================================================================================================    
% INPUT PROCESSING
% ===================================================================================================================


% DEFAULT ARGUMENTS
    function defInputs(paramSet)
        if ~paramSet(1) ; DIMS_K = find(size(Signal)~=1,1,'last') ; end
        if ~paramSet(2) ; FUNC = repmat('exp',[length(DIMS_K) 1]) ; end
        if ~paramSet(3) ; R0 = 1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)) ; end
        if ~paramSet(4) ; CRIT_THRS = 1 ; end
        if ~paramSet(5) ; FIT = 'TLS' ; end
        if ~paramSet(6) ; DECIM = ones(1,ndims(Signal)) ; end
        if ~paramSet(7) ; SHIFTS = eye(length(DIMS_K)) ; end
        if ~paramSet(8) ; CHOICE = 'auto' ; end
        if ~paramSet(9) ; STABILDIAG = false ; end
        if ~paramSet(10) ; MAC = false ; end
        if ~paramSet(11) ; DEBUG = false ; end
        if ~paramSet(12) ; CRITERION = 'MDL' ; end
        if ~paramSet(13) ; SOLVER = 'eig' ; end
        if ~paramSet(14) ; COMPUTE_U = false ; end
        if ~paramSet(15) ; M_L = 2/3 ; end
        if ~paramSet(16) ; W0 = [] ; end
        if ~paramSet(17) ; COMPUTE_dK = false ; end
        if ~paramSet(18) ; SIGNAL_MODEL = false ; end
        if ~paramSet(19) ; COMPUTE_dU = false ; end
        if ~paramSet(20) ; Kstab = [] ; end
        if ~paramSet(21) ; Ustab = [] ; end
    end


% PROCESS INPUTS
    function parseInputs
        % Initialize
            nargin = length(varargin) ;
            paramSet = false(100) ;
            defInputs(paramSet) ;
        % Processing to do
            TODO = [] ; % Default, DO everything
        % Is the first argument a structure ?
            if isstruct(varargin{1}) && strcmp(varargin{1}.Title,'ESPRIT_Results')
                % Initialize the output
                    OUT = varargin{1} ;
                % Import the structure
                    for fname = fieldnames(OUT)'
                        try eval([fname{1},'=OUT.',fname{1},';']) ; end
                    end
                % Is there a second argument to limit the processing ?
                    if nargin>2
                        TODO = varargin{2} ;
                    end
                % Skeep varargin processing
                    return ;
            end
        % Is DIMS_K the first argument ?
            if mod(nargin,2)==1 
                if isnumeric(varargin{1})
                    DIMS_K = varargin{1} ;
                    paramSet(1) = true ;
                else
                    errorInput('Wrong second argument : should be DIMS_K or a string') ;
                end
                if nargin>1
                    varargin = varargin(2:end) ;
                end
            end
        % Treat following arguments
            if nargin>1
                for i = 1:2:length(varargin)-1
                    Name = varargin{i} ;
                    Value = varargin{i+1} ;
                    if isempty(Value) ; continue ; end
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
                        case 'CRIT_THRS'
                            CRIT_THRS = Value ;
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
                        case 'MAC'
                            MAC = Value ;
                            paramSet(10) = true ;
                        case 'DEBUG'
                            DEBUG = Value ;
                            paramSet(11) = true ;
                        case 'CRITERION'
                            CRITERION = Value ;
                            paramSet(12) = true ;
                        case 'SOLVER'
                            SOLVER = Value ;
                            paramSet(13) = true ;
                        case 'COMPUTE_U'
                            COMPUTE_U = Value ;
                            paramSet(14) = true ;
                        case 'M/L'
                            M_L = Value ;
                            paramSet(15) = true ;
                        case 'W0'
                            W0 = Value ;
                            paramSet(16) = true ;
                        case 'COMPUTE_DK'
                            COMPUTE_dK = Value ;
                            paramSet(17) = true ;
                        case 'SIGNAL_MODEL'
                            SIGNAL_MODEL = Value ;
                            paramSet(18) = true ;
                        case 'COMPUTE_DU'
                            COMPUTE_dU = Value ;
                            paramSet(19) = true ;
                        case 'KSTAB'
                            Kstab = Value ;
                            paramSet(20) = true ;
                        case 'USTAB'
                            Ustab = Value ;
                            paramSet(21) = true ;
                        otherwise
                            %errorInput(['Wrong argument name in n°',num2str(i),'.'])
                            errorInput([Name,' (n°',num2str(i),').'])
                    end
                end
            end
        % DEFAULT VALUES
            defInputs(paramSet) ;
    end


% PROMPT AN ERROR ON WRONG INPUT ARGUMENTS
    function errorInput(info) 
        msg = ['Incorrect input argument : '] ;
        msg = [msg,info] ;
        error([msg,char(10)])
    end


    

% ===================================================================================================================    
% STABILIZATION DIAGRAM
% ===================================================================================================================

% PLOT THE STABILIZATION DIAGRAM
    function stabilizationDiagram()
        if length(DIMS_K)>1 ; warning('Stabilization Diagram is available for 1D-ESPRIT only.') ; return ; end
        % Are old results available or compute the stab. diag. ?
            if onlyStabDiag
                K = Kstab(end,:) ;
                R = size(Kstab,1) ;
            else
                % Compute All the Poles for All Signal Orders
                    Kstab = zeros(length(R0),max(R0))*NaN*(1+1i) ;
                    for r = 1:length(R0)
                        extractPoles(r) ; 
                        Kstab(r,1:R0(r)) = K ;
                    end
                    K = Kstab(R,1:R0(R)) ;
                % Sort the poles
                    Kstab = sort(Kstab,2,'ascend') ; % Simply sort by value
            end
        % Initialize the stabilization diagram
            stab = InitStabDiag() ;
        % SORT WITH MAC VALUES
            % Initialize the branches
                stab.branches = plot3(stab.axPoles,...
                                        NaN*ones(length(R0),max(R0)),...
                                        NaN*ones(length(R0),max(R0)),...
                                        NaN*ones(length(R0),max(R0)),...
                                        '-k','linewidth',.5) ;
                uistack(stab.branches,'bottom') ;
                set(stab.branches,'color',[1 1 1]*.5) ;
            % Sort if needed
                if isempty(Ustab) && MAC % Sort with MACs
                    stab = SortStabDiag(stab) ;
                    V = [] ;
                end
        % Backup the current poles choice
            K = Kstab(R,1:R0(R)) ;
        % Init interactions
            stab = InitModeShape(stab) ;
            stab = InitStabDiagInteract(stab) ;
        % Wait for the figure to be closed
            stab;
            if ~onlyStabDiag ; uiwait(stab.fig) ; end
            drawnow ;
    end


% FOLLOW MODES WITH MAC VALUES
    function stab = SortStabDiag(stab)    
        % Initialize
            %Ustab = ones(prod(Lp),length(R0),max(R0))*NaN*(1+1i) ; % dims : [Point order pole]
            Ustab = ones(length(indP),length(R0),max(R0))*NaN*(1+1i) ; % dims : [Point order pole]
        % Process Data
            for r = 1:length(R0)
                % Compute the mode
                    K = Kstab(r,1:R0(r)) ;
                    buildVandermonde(Lk) ;
                    Ustab(:,r,1:R0(r)) = (V\Signal(indP,:).').' ;
                % Not the first mode ?
                    if r>1
                        % Compute the MAC matrix
                            U1 = squeeze(Ustab(:,r-1,1:R0(r-1))) ;
                            U2 = squeeze(Ustab(:,r,1:R0(r))) ;
                            mac = abs(U2'*U1)./sqrt(sum(abs(U2).^2,1)'*sum(abs(U1).^2,1)) ;
                        % Sort the parameters
                            indices = NaN*ones(1,R0(r-1)) ;
                            % Loop
                                for rr = 1:R0(r-1)
                                    % Get the max. Value
                                        maxMAC = max(mac(:)) ;
                                        if maxMAC<MAC ; break ; end
                                        [rmax,cmax] = find(mac==maxMAC) ;
                                    % Set the corresponding mac values to NaN
                                        mac(:,cmax(1)) = nan ;
                                        mac(rmax(1),:) = nan ;
                                    % Save to the indices vector
                                        indices(cmax(1))=rmax(1) ;
                                    % Break if all modes have been classified
                                        if ~any(~isnan(mac(:))) ; break ; end
                                end
                            % Remaining non-sorted indices
                                indices = indices(~isnan(indices)) ;
                                indices = [indices setdiff(1:R0(r),indices)] ;
                            % Rearrange
                                Kstab(r,1:R0(r)) = Kstab(r,indices) ;
                                Ustab(:,r,1:R0(r)) = Ustab(:,r,indices) ;
                    end
                % Update Plot
                    for br = 1:R0(r)
                        stab.branches(br).XData(1:r) = abs(real(Kstab(1:r,br))) ;
                        stab.branches(br).YData(1:r) = R0(1:r) ;
                        stab.branches(br).ZData(1:r) = abs(imag(Kstab(1:r,br))./real(Kstab(1:r,br))) ;
                    end
                % Update the sorted scatter plot
                    stab.scatCriterion.XData = abs(real(Kstab(:))) ;
                    stab.scatCriterion.ZData = abs(imag(Kstab(:))./real(Kstab(:))) ;
                % Draw
                    drawnow ;
            end
    end
    


% STABILIZATION DIAGRAM INITIALIZATION
    function stab = InitStabDiag()
        % Figure init
            relSize = .85 ;
            stab.fig = figure(...
                            'NumberTitle','off'...
                            ,'Name','ESPRIT : STABILIZATION DIAGRAM. Click to select an other signal order. Right Click to quit.'...
                            ,'tag','ESPRIT.StabDiag' ...
                            ) ;
            %if DEBUG ; stab.fig.WindowStyle = 'docked' ; end
            stab.fig.Position(1:2) = stab.fig.Position(1:2) + stab.fig.Position(3:4)*(1-relSize)/2 ;
            stab.fig.Position(3:4) = stab.fig.Position(3:4)*relSize ;
        % Axes for the Poles
            stab.axPoles = axes('outerposition',[.2 0 .8 1]) ;
                stab.scatCriterion = scatter3(stab.axPoles,...
                        abs(real(Kstab(:))),...
                        repmat(R0(:),[size(Kstab,2),1]),...
                        abs(imag(Kstab(:))./real(Kstab(:))),...
                        200,...
                        'k',...
                        '.',...
                        'tag','scatCriterion') ;
                %plot3(abs(real(Kstab)),repmat(R0(:),[1 max(R0)]),abs(imag(Kstab)./real(Kstab)),'.k','markersize',10,'linewidth',.5)
                stab.axPoles.ZScale = 'log' ;
                stab.axPoles.SortMethod = 'childorder' ;
                % Format
                    stab.axPoles.YTickLabel = []; 
                    stab.axPoles.YMinorTick = 'on'; 
                    %stab.axPoles.YMinorGrid = 'on'; 
                    grid(stab.axPoles,'on') ;
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
                % Limits
                    %set([stab.axMean stab.axPoles],'xlim',[1/prod(Lk)*2*pi/2 pi]) ;
                    w = abs(real(Kstab(:))) ;
                    set([stab.axMean stab.axPoles],'xlim',[min(w(w~=0)) max(w)].*[0.9 1.1]) ;
                set([stab.axPoles],'ylim',[min(R0) max(R0)]) ;
                global hlink , hlink = linkprop([stab.axMean stab.axPoles],'position') ;
                linkaxes([stab.axMean stab.axPoles],'x') ;
                uistack(stab.axMean,'bottom') ;
        % Axes for the signal order selection Criterion(s)
            stab.axCrit = axes('outerposition',[0 0 .2 1]) ;
                % Plot ALL Criterions
                    %plot(R0,(CRIT./max(CRIT)),'.-','markersize',20,'linewidth',1) ;
                    plot(R0,(MDL./max(MDL)),'.-','markersize',20,'linewidth',1) ;
                    plot(R0,(ESTER./max(ESTER)),'.-','markersize',20,'linewidth',1) ;
                    plot(R0,(SAMOS./max(SAMOS)),'.-','markersize',20,'linewidth',1) ;
                % Format the axes
                    set(gca,'view',[-90 90]) ;
                    box on
                    grid on
                    axis tight
                    stab.axCrit.XLim = [min(R0) max(R0)] ;
                % Cursors
                    stab.plOrderLine = plot(stab.axCrit,R0(R)*[1 1],stab.axCrit.YLim,'-.k','linewidth',1) ;
                    stab.plRLine = plot(stab.axCrit,R0(R)*[1 1],stab.axCrit.YLim,'-.r','linewidth',1) ;
                % Disable rotate3d
                    hBehavior = hggetbehavior(stab.axCrit,'Rotate3d');
                    hBehavior.Enable = false ;
                % Legend
                    legend({'MDL','ESTER','SAMOS'},'edgecolor','k','location','northwest') ;
                % Re-align with the poles axes
                    stab.axCrit.Position([2,4]) = stab.axPoles.Position([2,4]) ;
        % Poles Chosen at the end
            stab.plRPoles = plot3(stab.axPoles,abs(real(K)),R0(R)*ones(1,R0(R)),abs(imag(K)./real(K)),'or','markersize',15,'linewidth',1.5) ;
        % Poles at the mouse position
            stab.plOrderPoles = plot3(stab.axPoles,abs(real(K)),R0(R)*ones(1,R0(R)),abs(imag(K)./real(K)),'.m','markersize',25) ;
        % DRAW
            drawnow ;
    end


% MODE SHAPE PLOT
    function stab = InitModeShape(stab) 
        % Prepare the axes
            stab.axShape = axes() ;
            stab.axShape.Position(1:2) = stab.axPoles.Position(1:2) ;
            stab.axShape.Position(3:4) = [.2 .2] ;
            stab.axShape.XTick = [] ;
            stab.axShape.YTick = [] ;
            stab.axShape.XLim = [-1 1] ;
            stab.axShape.YLim = [-1 1] ;
            axis(stab.axShape,'equal') ;
            box(stab.axShape,'on') ;
        % Prepare the plot
            stab.plotShape = plot(stab.axShape,NaN,NaN,'+b','linewidth',.5,'markersize',8) ;
            stab.markShape = plot(stab.axPoles,NaN,NaN,'+b','linewidth',2,'markersize',15) ;
    end


% INITIALIZE THE SAB DIAG INTERACTION
    function stab = InitStabDiagInteract(stab)
        % Buttons for the view
            margin = 0.003 ;
            btnWidth = 0.08 ;
            btnHeight = 0.03 ;
            % Frequency/damping switch
                stab.btnSwitch = uicontrol(stab.fig,'style','pushbutton') ;
                stab.btnSwitch.String = 'Frequency' ;
                stab.btnSwitch.TooltipString = 'Switch Representation Mode' ;
                stab.btnSwitch.Units = 'normalized' ;
                stab.btnSwitch.Position = [stab.axPoles.Position(1:2)+stab.axPoles.Position(3)*[1 0]+[-btnWidth-margin btnHeight+margin] [btnWidth btnHeight]] ;
                stab.btnSwitch.Callback = @(src,evt)btnSwitchCallback(stab) ;
            % Listbox choose selection mode
                stab.popupSelect = uicontrol(stab.fig,'style','popupmenu') ;
                stab.popupSelect.String = {'Full','One'} ;
                stab.popupSelect.TooltipString = 'Pole Selection' ;
                stab.popupSelect.Units = 'normalized' ;
                stab.popupSelect.Position = [stab.axPoles.Position(1:2)+stab.axPoles.Position(3)*[1 0]+[-2*btnWidth-2*margin btnHeight+margin] [btnWidth btnHeight]] ;
                stab.popupSelect.Callback = @(src,evt)popupSelectCallback(stab) ;
            % Listbox Criterion plot
                stab.popupCriterion = uicontrol(stab.fig,'style','popupmenu') ;
                stab.popupCriterion.String = {'None','Complexity','MAC','sigma(f)','sigma(xi)','sigma(phi)'} ;
                stab.popupCriterion.TooltipString = 'Show a criterion' ;
                stab.popupCriterion.Units = 'normalized' ;
                stab.popupCriterion.Position = [stab.axPoles.Position(1:2)+stab.axPoles.Position(3)*[1 0]+[-3*btnWidth-3*margin btnHeight+margin] [btnWidth btnHeight]] ;
                stab.popupCriterion.Callback = @(src,evt)popupCriterionCallback(stab) ;
            % Slider to tune the criterion
                stab.sliderCriterion = uicontrol(stab.fig,'style','slider','tag','sliderCriterion') ;
                stab.sliderCriterion.TooltipString = 'Tune the criterion' ;
                stab.sliderCriterion.Units = 'normalized' ;
                stab.sliderCriterion.Position = [stab.axPoles.Position(1:2)+stab.axPoles.Position(3)*[1 0]+[-5*btnWidth-4*margin btnHeight+margin] [2*btnWidth btnHeight]] ;
                addlistener(stab.sliderCriterion, 'ContinuousValueChange', @(src,evt) sliderCriterionCallback(stab)) ;
                %stab.sliderCriterion.Callback = @(src,evt)sliderCriterionCallback(stab) ;
        % Figure Callbacks setting
            stab.fig.WindowButtonMotionFcn = @(src,evt)changeStabDiagOrder(Kstab,stab,'move') ;
            stab.fig.WindowButtonDownFcn = @(src,evt)changeStabDiagOrder(Kstab,stab,'click') ;
    end

% CHANGE THE STABILIZATION DIAGRAM ORDER WITH MOUSE POSITION
    function changeStabDiagOrder(Kstab,stab,event)
        % Get the order given by the mouse position
            mouseK = stab.axPoles.CurrentPoint(1,1) ;
            mouseOrder = stab.axCrit.CurrentPoint(1,1) ;
            selectOrder = min(max(R0(1),round(mouseOrder)),R0(end)) ;
        % Get the pole
            [~,order] = min(abs(R0-selectOrder)) ;
            kk = abs(real(Kstab(order,:))) ;
            [~,numPole] = min(abs(kk-mouseK)) ;
        % Do something
        switch event % MOUSE MOVED OR CLICKED
            case 'move' % ONLY DISPLAY
                % Process the mode shape
                    if MAC
                        uu = Ustab(:,order,numPole) ;
                        uu = uu./max(abs(uu(:))) ;
                        % Plot the mode shape
                            stab.plotShape.XData = real(uu) ;
                            stab.plotShape.YData = imag(uu) ;
                            stab.markShape.XData = abs(real(Kstab(order,numPole))) ;
                            stab.markShape.YData = R0(order) ;
                            stab.markShape.ZData = abs(imag(Kstab(order,numPole))./real(Kstab(order,numPole))) ;
                    end
                % Order cursor
                    stab.plOrderLine.XData = R0(order)*[1 1] ;
                    stab.plOrderPoles.XData = abs(real(Kstab(order,:))) ;
                    stab.plOrderPoles.YData = R0(order)*ones(1,max(R0)) ;
                    stab.plOrderPoles.ZData = abs(imag(Kstab(order,:))./real(Kstab(order,:))) ;
            case 'click' % CHANGE THE CHOOSEN POLES
                switch stab.fig.SelectionType % Mode of selection
                    case 'normal' % ADD THE CORRRESPONDING POLES
                        % Depending on the selection mode...
                            switch stab.popupSelect.String{stab.popupSelect.Value}
                                case 'Full'
                                    R = order ;
                                    K = Kstab(R,1:R0(R)) ;
                                    stab.plRPoles.XData = abs(real(K)) ;
                                    stab.plRPoles.YData = R0(order)*ones(1,R0(R)) ;
                                    stab.plRPoles.ZData = abs(imag(K)./real(K)) ;
                                case 'One'
                                    % Has the pole already been selected ?
                                        ind = find(abs(abs(real(K))-abs(real(Kstab(order,numPole))))<eps) ;
                                        switch isempty(ind)
                                            case true % The pole is new: add it
                                                if DEBUG ; disp(['Add K: ',num2str(Kstab(order,numPole),3)]) ; end
                                                if ~isreal(Signal) % The signal is complex
                                                    K(end+1) = Kstab(order,numPole) ;
                                                    stab.plRPoles.YData(end+1) = R0(order) ;
                                                else % The signal is real, add the conjugate pole
                                                    K(end+(1:2)) = [Kstab(order,numPole),-conj(Kstab(order,numPole))] ;
                                                    stab.plRPoles.YData(end+(1:2)) = R0(order)*[1 1] ;
                                                end 
                                            case false % The pole has to be removed
                                                if DEBUG ; disp(['Rem K: ',num2str(Kstab(order,numPole),3)]) ; end
                                                K(ind) = [] ;
                                                stab.plRPoles.YData(ind) = [] ;
                                        end
                                    % Update the markers
                                        stab.plRPoles.XData = abs(real(K)) ;
                                        stab.plRPoles.ZData = abs(imag(K)./real(K)) ;
                            end
                        % Update markers
                            stab.plRLine.XData = R0(order)*[1 1] ;
                        % DISPLAY INFOS
                            if DEBUG ; disp(['K: ',mat2str(K,2)]) ; end
                    case 'open' % NO ACTION
                    case 'alt' % RIGHT-CLICK, CLOSE THE STAB. DIAG.
                        close(stab.fig) ;
                end
        end
    end

% CHANGE PLOT CRITERION
    function popupCriterionCallback(stab)
        % Delete any preovious criterion
            %delete(findobj(stab.axPoles,'tag','scatCriterion')) ;
            delete(findobj(stab.fig,'type','colorbar')) ;
        % Compute the criterion
            crit = NaN*ones(size(Kstab)) ;
            maxCrit = 9 ;
            switch stab.popupCriterion.String{stab.popupCriterion.Value}
                case 'None' % No Criterion
                    stab.CData = 'k' ;
                    set(stab.branches,'visible','on') ;
                case 'Complexity' % Evaluate the complexity
                    for or = 1:size(Kstab,1) 
                        for rr = 1:size(Kstab,2) ;
                            if isnan(Kstab(or,rr)) ; continue ; end
                            uu = Ustab(:,or,rr) ;
                            uu = uu./abs(uu) ;
                            S = [real(uu(:)) imag(uu(:))] ;
                            g = sqrt(real(eig(S'*S))) ;
                            crit(or,rr) = min(g)/max(g) ;
                        end
                    end
                    maxCrit = 3 ;
                case 'MAC' % MAC Value between successive modes
                    for or = 2:size(Kstab,1) 
                        for rr = 1:size(Kstab,2) ;
                            if isnan(Kstab(or,rr)) ; continue ; end
                            if isnan(Kstab(or-1,rr)) ; continue ; end
                            mac = Ustab(:,or,rr)'*Ustab(:,or-1,rr)/norm(Ustab(:,or,rr))/norm(Ustab(:,or-1,rr)) ;
                            crit(or,rr) = 1-abs(mac) ;
                        end
                    end
                    maxCrit = 3 ;
                case 'sigma(phi)' % MAC Value between successive modes
                    for or = 2:size(Kstab,1) 
                        for rr = 1:size(Kstab,2) ;
                            if isnan(Kstab(or,rr)) ; continue ; end
                            if isnan(Kstab(or-1,rr)) ; continue ; end
                            crit(or,rr) = 1/2*norm(Ustab(:,or,rr)-Ustab(:,or-1,rr))/max(norm(Ustab(:,or,rr)),norm(Ustab(:,or-1,rr))) ;
                        end
                    end
                case 'sigma(f)' % MAC Value between successive modes
                    freqs = abs(real(Kstab)) ;
                    for or = 2:size(Kstab,1) 
                        for rr = 1:size(Kstab,2) ;
                            if isnan(Kstab(or,rr)) ; continue ; end
                            if isnan(Kstab(or-1,rr)) ; continue ; end
                            crit(or,rr) = abs(freqs(or,rr)-freqs(or-1,rr))/max(freqs(or,rr),freqs(or-1,rr)) ;
                        end
                    end
                case 'sigma(xi)' % MAC Value between successive modes
                    xi = abs(imag(Kstab)./real(Kstab)) ;
                    for or = 2:size(Kstab,1) 
                        for rr = 1:size(Kstab,2) ;
                            if isnan(Kstab(or,rr)) ; continue ; end
                            if isnan(Kstab(or-1,rr)) ; continue ; end
                            crit(or,rr) = abs(xi(or,rr)-xi(or-1,rr))/max(xi(or,rr),xi(or-1,rr)) ;
                        end
                    end
            end
        % Scatter plot
            if ~strcmp('None',stab.popupCriterion.String{stab.popupCriterion.Value})
                % Hide the branches
                    set(stab.branches,'visible','off') ;
                % Scatter
                    stab.scatCriterion.CData = -log10(crit(:)) ;
                % Colorbar
                    colorbar(stab.axPoles) ;
                    colormap(linspecer(1000)*.85) ;
                % Color Range
                    caxis auto
                    caxis(stab.axPoles,[min(caxis(stab.axPoles)) min(max(caxis(stab.axPoles)),maxCrit)])
                % Slider range
                    crange = caxis(stab.axPoles) ;
                    sli = findobj(stab.fig,'tag','sliderCriterion') ;
                    sli.Min = crange(1) ;
                    sli.Max = crange(2) ;
                    sli.Value = sli.Min ;
            end
    end

% CHANGE PLOT CRITERION
    function sliderCriterionCallback(stab)
        scat = findobj(stab.axPoles,'tag','scatCriterion') ;
        if isempty(scat) ; return ; end
        cri = scat.CData ;
        fr = abs(real(Kstab(:))) ;
        fr(cri<=stab.sliderCriterion.Value) = NaN ;
        set(scat,'XData',fr) ;
    end

% CHANGE SELECTION MODE
    function popupSelectCallback(stab)
        % Initialize the choice of K
            switch stab.popupSelect.String{stab.popupSelect.Value}
                case 'Full' % Select the full order
                    K = Kstab(R,1:R0(R)) ;
                    stab.plRPoles.YData = R0(R)*ones(1,R0(R)) ;
                case 'One' % Select poles one by one
            end
        % Initialize cursors
            stab.plRLine.XData = R0(R)*[1 1] ;
            stab.plRPoles.XData = abs(real(K)) ;
            stab.plRPoles.ZData = abs(imag(K)./real(K)) ;
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

