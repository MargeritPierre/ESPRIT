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
%           string : 'ESTER' or 'MDL' ('MDL')
%       CRIT_THRS : signal order criterion threshold 
%           1x1 double (1)
%       COMPUTE_U : compute amplitudes 
%           boolean (false)
%       SIGNAL_MODEL : compute the signal model 
%           boolean (false)
%       COMPUTE_dK : compute uncertainties 
%           boolean (false)
%       DIMS_K : dimensions of the wavevectors 
%           1xD integer (ndims(Signal))
%       FUNC : type of function searched 
%           Dx3 char : 'exp' or 'cos' (repmat('exp',[length(DIMS_K) 1]))
%       R0 : signal order candidates 
%           1xnR integer : (1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)))
%       M/L : Spatial Smoothing ratio (gives the HbH matrix shape)
%           1xD double : 0<M_d/L<1 (2/3)
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

    % Input Initialization
        parseInputs ; 
        varargin ;
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
        
   % Output initialization
        OUT = [] ;

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
        K = zeros([max(R0) length(DIMS_K)]) ; % Wavevectors Matrix
        dK = zeros([max(R0) length(DIMS_K)]) ; % Estimated Variance
        V = [] ; % Steering Matrix
        U = zeros([max(R0) length(DIMS_P)+1]) ; % Amplitudes
        SignalModel = [] ; % Reconstructed Signal Model
        Hss = [] ;
        
        
        
        
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
                Hss = buildHss(Signal) ;
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
    % KEEP the needed E.V only
        lambda = lambda(1:max(R0)) ;
        W = W(:,1:max(R0)) ;
    % GET THE LEFT EIGENVECTORS
        if transpCss
            W = Hss*W ;
            W = bsxfun(@(w,l)w./l, W, sqrt(lambda(:)).') ;
        end
        
        
% SIGNAL ORGER CRITERION
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Signal Order Criterion : ') ; lastTime = tic ; end
    if length(R0)==1 % THE SIGNAL ORDER HAS BEEN GIVEN
        R = 1 ;
    else % SIGNAL ORDER CHOICE
        % COMPUTE THE CRITERION
            switch CRITERION
                case 'ESTER'
                    % Compute all errors
                        ESTER = zeros(size(SHIFTS,1),length(R0)) ;
                        for r = 1:length(R0)
                            for s = 1:size(SHIFTS,1)
                                ESTER(s,r) = ester(r,s) ;
                            end
                        end
                    % Compute the mean over shifts
                        %CRIT = prod(ESTER,1) ;
                        CRIT = max(ESTER,[],1) ;
                        %CRIT = mean(ESTER,1) ;
                    % Save the criterion
                        OUT.ESTER = ESTER ;
                case 'MDL'
                    MDL = (lambda(1:end-1)-lambda(2:end))./lambda(1:end-1) ;
                    CRIT = [MDL ; 0] ;
                    CRIT = CRIT(R0) ;
                    % Save the criterion
                        OUT.MDL = MDL ;
            end
            % CHOOSE THE SIGNAL ORDER
                R = find(CRIT>=max(CRIT)*CRIT_THRS,1,'last') ;
                OUT.CRIT = CRIT ;
    end
    
        
% STABILIZATION DIAGRAM
    if STABILDIAG
        MAC ;
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Stabilization Diagram : ') ; lastTime = tic ; end
        stabilizationDiagram() ;
    end
        
    
% EXTRACT THE WAVEVECTORS
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Pole Estimation : ') ; lastTime = tic ; end
    [T,beta] = extractPoles(R) ;
        
    
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
    if COMPUTE_dK 
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Uncertainty Estimation : ') ; lastTime = tic ; end
        computeUncertainties ;
    end
        
    
% OUTPUT
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Output : ') ; lastTime = tic ; end
    % Input arguments
        OUT.varargin = varargin ;
    % Results
        OUT.K = K ;
        OUT.dK = dK ;
        OUT.V = V ;
        OUT.U = U ;
        if SIGNAL_MODEL
            if isempty(DIMS_P)
                OUT.SignalModel = reshape(SignalModel,Lk) ;
            else
                OUT.SignalModel = permute(reshape(SignalModel,[Lp Lk]),permDims) ;
            end
        end
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
                        if p==1 && DEBUG ; disp(['       estim: ',num2str(toc(t)*length(indP)),' secs']) ; end
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
    function [F,Wup,Wdwn] = computeF(r,s)
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
        % Shift invariance evaluation
            switch FIT
                case 'LS'
                    F = Wup\Wdwn ;
                case 'TLS'
                    % Uncorrelated data
                        kTLS = floor(size(Wp,1)/2) ;
                        Wup = Wp(2*(0:kTLS-1)+1,:) ;
                        Wdwn = Wp(2*(0:kTLS-1)+2,:) ;
                    % TLS
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
    function [T,Beta] = extractPoles(r)
        % Evaluation of all shift invariances
            F = zeros(R0(r),R0(r),size(SHIFTS,1)) ;
            for s = 1:size(SHIFTS,1)
                F(:,:,s) = computeF(r,s) ;
            end
        % Evaluation of PHI
            if size(SHIFTS,1)==1 % One shift, simple diagonalization
                [T,PHI] = eig(F) ;
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
        % Poles
            Z = exp(1i.*K) ;
            Z = repmat(reshape(Z,[1 size(Z)]),[prod(L) 1]) ;
        % Vandermonde
            indV = repmat(indV,[1 1 size(Z,3)]) ;
            V = reshape(prod(Z.^indV,2),[prod(L) size(Z,3)]) ;
    end


% AMPLITUDES U RETRIEVAL
    function computeU % Here, decimation factors are no longer taken into account
        % Vandermonde Matrix
            buildVandermonde(Lk) ;
        % Shift for the conditionning
%             v0 = max(abs(V),[],1) ;
%             V = V*diag(1./v0) ;
        % Signal reshaping
            S = Signal.' ;
        % Amplitude estimation
            A = V\S ;
        % Vandermonde Shift compensation
%             A = diag(1./v0)*A ;
        % Reshaping
            U = reshape(A.',[Lp , size(K,2)]) ;
    end

% SIGNAl MODEL CONSTRUCTION
    function computeSignalModel()
        if isempty(V) ; computeU() ; end
        SignalModel = (V*A).' ;
    end


% UNCERTAINTY ESTIMATION dK
    function computeUncertainties
        % Options (hard-coded for now)
            lin_method = 'conv' ; % Linearization method for the variance: 'none', 'kron' or 'conv'
            covar_estim = 'uniform' ; % estimation of the perturbation covariance: 'uniform', 'diagonal' or 'full'
        % Init
            dK = zeros(size(K)) ;
        % Signal perturbation
            if isempty(SignalModel) ; computeSignalModel() ; end
            dS = Signal-SignalModel ;
        % Perturbation covariance
            if ~strcmp(lin_method,'none')
                switch covar_estim
                    case 'uniform'
                        var_dS = var(dS(:),0);%*speye(prod(Lk)*length(indP)) ;
                    case 'diagonal'
                        var_dS = diag(var(dS,0,1)) ;
                    case 'full'
                        dSzm = (dS-repmat(mean(dS,1),[size(dS,1) 1])).' ;
                        var_dS = dSzm*dSzm'/size(dSzm,2) ;
                end
            end
        % Selection matrices if needed
            switch lin_method
                case 'none'
                    dH = buildHss(dS) ;
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
                case 'conv'
                    % Jdelta matrix (decimation selection matrix)
                        I = speye(prod((Kk-1).*DECIM_K+1)) ;
                        Jdelta = I(indHbH{1}(:,1),:) ;
            end
        % Pre-Computations
            shifts = eye(length(DIMS_K)) ; % /!\ SHIFTS IS DEFAULT HERE !
            % Partial Vandermonde matrices
                if isempty(V) ; buildVandermonde(Lk) ; end
                P = V(indHbH{1}(:,1),:) ;
                Q = V(indHbH{1}(1,:),:) ;
            % Complete right-Vandermonde Matrix
                QA = zeros(prod(Mk)*length(indP),R0(R)) ;
                for p = 1:length(indP) 
                    QA((1:prod(Mk))+(p-1)*prod(Mk),:) = Q*diag(A(:,indP(p))) ;
                end
            % Right vector of the bilinear form
                if 0
                    x = (QA\eye(size(QA,1)))' ;
                else
                    if isempty(Hss) ; Hss = buildHss(SignalModel) ; end
                    x = conj(bsxfun(@(x,c)x./c,Hss'*W(:,1:R0(R)),lambda(1:R0(R)).')*T) ;
                end
        % Loop over the shifts
            for s = 1:size(K,1)
                [~,~,~,Jup,Jdwn] = selectMatrices(shifts(s,:)) ;
                % pre-compute
                    if 0
                        vn = (Jup*P)\speye(size(Jdwn,1)) ;
                    else
                        vn = T\((Jup*W(:,1:R0(R)))\speye(size(Jdwn,1))) ;
                    end
                % Loop over the R0(R) components
                    for r = 1:size(K,2)
                        % The pole
                            arn = P(2,r) ;
                        % Left vector of the bilinear form
                            vrn = (vn(r,:)*(Jdwn-arn*Jup)*Jdelta)' ;
                        % Linearization method
                            switch lin_method
                                case 'none' % UNCERTAINTY (DOES NOT WORK WELL...)
                                    dK(s,r) = abs(vrn'*dH*x(:,r)/arn/prod(DECIM_K))^2 ;
                                case 'conv' % LINEAR / BY CONVOLUTION (USE OF THE HANKEL SHAPE OF Hss)
                                    VRN = reshape(vrn,[1 (Kk-1).*DECIM_K+1]) ;
                                    XR = reshape(x(:,r),[length(indP) Mk]) ;
                                    ZN = convn(VRN,XR) ;
                                    zn = ZN(:) ;
                                case 'kron' % LINEAR / BY VECTORIZATION (USES vec(A*X*B) = kron(B.',A)*vec(X) ) 
                                    zn = (kron((x(:,r)),vrn)'*M).' ;
                            end
                        % Save
                            dK(s,r) = zn'*var_dS*zn/prod(DECIM_K)^2 ;
                    end
            end
    end


    

% ===================================================================================================================    
% FUNCTIONS FOR UTILS
% ===================================================================================================================

% PROCESS INPUTS
    function parseInputs
        nargin = length(varargin) ;
        paramSet = false(100) ;
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
                        case 'CRITERION'
                            CRITERION = Value ;
                            paramSet(12) = true ;
                        case 'CRIT_THRS'
                            CRIT_THRS = Value ;
                            paramSet(4) = true ;
                        case 'COMPUTE_U'
                            COMPUTE_U = Value ;
                            paramSet(14) = true ;
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
                        case 'SOLVER'
                            SOLVER = Value ;
                            paramSet(13) = true ;
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
                        otherwise
                            %errorInput(['Wrong argument name in n°',num2str(i),'.'])
                            errorInput([Name,' (n°',num2str(i),').'])
                    end
                end
            end
        % DEFAULT VALUES
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
    end


% PROMPT AN ERROR ON WRONG INPUT ARGUMENTS
    function errorInput(info) 
        msg = ['Incorrect input argument : '] ;
        msg = [msg,info] ;
        error([msg,char(10)])
    end


    

% ===================================================================================================================    
% GRAPHICAL FUNCTIONS
% ===================================================================================================================

% PLOT THE STABILIZATION DIAGRAM
    function stabilizationDiagram()
        if length(DIMS_K)>1 ; warning('Stabilization Diagram is available for 1D-ESPRIT only.') ; return ; end
        % Compute All the Poles for All Signal Orders
            Kstab = zeros(length(R0),max(R0))*NaN*(1+1i) ;
            for r = 1:length(R0)
                extractPoles(r) ; 
                Kstab(r,1:R0(r)) = K ;
            end
        % Sort the poles
            if ~MAC % Simply sort by value
                Kstab = sort(Kstab,2) ;
            else % Sort with MACs
                % Compute all Modes for all Mode Orders
                    Ustab = zeros(prod(Lp),length(R0),max(R0))*NaN*(1+1i) ; % dims : [Point order pole]
                    wtbr = waitbar(0,'Computing Modes...') ;
                    for r = 1:length(R0)
                        K = Kstab(r,1:R0(r)) ;
                        computeU ;
                        Ustab(:,r,1:R0(r)) = reshape(U,[prod(Lp),R0(r)]) ;
                        wtbr = waitbar(r/length(R0),wtbr) ;
                    end
                    delete(wtbr) ; drawnow ;
                % Compute MAC Values
                    K_MAClinks = [] ;
                    OR_MAClinks = [] ;
                    for or = length(R0)-1:-1:1 ;
                        U1 = squeeze(Ustab(:,or+1,1:R0(or+1))) ;
                        U2 = squeeze(Ustab(:,or,1:R0(or))) ;
                        mac = abs(U2'*U1)./sqrt(sum(abs(U2).^2,1)'*sum(abs(U1).^2,1)) ;
                        for rr = 1:or
                            maxMAC = max(mac(:)) ;
                            if maxMAC<MAC ; break ; end
                            [rmax,cmax] = find(mac==maxMAC) ;
                            K_MAClinks = [K_MAClinks ; [Kstab(or+1,cmax(1)) Kstab(or,rmax(1))]] ;
                            OR_MAClinks = [OR_MAClinks ; [or+1 or]] ;
                            mac(:,cmax(1)) = nan ;
                            mac(rmax(1),:) = nan ;
                        end
                    end
            end
        % Backup the current poles choice
            K = Kstab(R,1:R) ;
        % Figure
            stab.fig = figure('NumberTitle','off','Name','ESPRIT : STABILIZATION DIAGRAM. Click to select an other signal order. Right Click to quit.') ; %figure ;
            % Axes for the signal order selection Criterions
                stab.axCrit = axes('outerposition',[0 0 .2 1]) ;
                    plot(R0,log10(CRIT),'.-','markersize',20,'linewidth',1) ;
                    set(gca,'view',[-90 90]) ;
                    box on
                    grid on
                    axis tight
                    stab.axCrit.XLim = [min(R0) max(R0)] ;
                    % Cursors
                        stab.plOrderLine = plot(stab.axCrit,R0(R)*[1 1],stab.axCrit.YLim,'-.k','linewidth',1) ;
                        stab.plRLine = plot(stab.axCrit,R0(R)*[1 1],stab.axCrit.YLim,'-.r','linewidth',1) ;
                    % Disable rotate3d
                        hBehavior = hggetbehavior(stab.axCrit, 'Rotate3d');
                        hBehavior.Enable = false ;
            % Axes for the Poles
                stab.axPoles = axes('outerposition',[.2 0 .8 1]) ;
                    if MAC ; plot3(abs(real(K_MAClinks')),OR_MAClinks',abs(imag(K_MAClinks')./real(K_MAClinks')),'-k','linewidth',.5) ; end
                    plot3(abs(real(Kstab)),repmat(R0(:),[1 max(R0)]),abs(imag(Kstab)./real(Kstab)),'+','markersize',10,'linewidth',.5)
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
                % Limits
                    %set([stab.axMean stab.axPoles],'xlim',[1/prod(Lk)*2*pi/2 pi]) ;
                    w = abs(real(Kstab(:))) ;
                    set([stab.axMean stab.axPoles],'xlim',[min(w(w~=0)) max(w)].*[0.9 1.1]) ;
                set([stab.axPoles],'ylim',[min(R0) max(R0)]) ;
                global hlink , hlink = linkprop([stab.axMean stab.axPoles],'position') ;
                linkaxes([stab.axMean stab.axPoles],'x') ;
                uistack(stab.axMean,'bottom') ;
            % Poles Chosen at the end
                stab.plRPoles = plot3(stab.axPoles,abs(real(Kstab(R,:))),R0(R)*ones(1,max(R0)),abs(imag(Kstab(R,:))./real(Kstab(R,:))),'or','markersize',15,'linewidth',1.5) ;
            % Poles at the mouse position
                stab.plOrderPoles = plot3(stab.axPoles,abs(real(Kstab(R,:))),R0(R)*ones(1,max(R0)),abs(imag(Kstab(R,:))./real(Kstab(R,:))),'.k','markersize',30) ;
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
                drawnow ;
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

