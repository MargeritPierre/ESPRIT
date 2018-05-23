function OUT = ESPRIT(Signal,varargin)
%
% ESPRIT Extract wavevectors from a signal
%
%   OUT = ESPRIT(Signal,varargin)
%
%   varargin parameters (and default values) : 
%       CHOICE : 'auto' or 'manual' ('auto')
%       SOLVER : 'eigs' or 'eig' ('eig') 
%       CRITERION : 'ESTER' or 'MDL' ('MDL')
%       CRIT_THRS : criterion threshold (1)
%       COMPUTE_U : boolean to compute amplitudes (true)
%       COMPUTE_dK : boolean to compute uncertainties (true)
%       DIMS_K : dimensions of the wavevectors (ndims(Signal))  
%       FUNC : type of function searched 'exp' or 'cos' (repmat('exp',[length(DIMS_K) 1]))
%       R0 : signal order candidates (1:floor(min(arrayfun(@(d)size(Signal,d),DIMS_K)/2)))
%       K0 : length(s) of the signal that gives the size of the autocovariance matrix ([] for an auto choice)
%       W0 : signal subspace of a previous iteration ([] for a random subspace initialization)
%       FIT : 'LS' or 'TLS' (TLS)
%       DECIM : decimation factors (ones(1,ndims(Signal)))
%           - along DIMS_K : /!\ NYQUIST
%           - along DIMS_P : reduce the size of Css. However, U is estimated at all points.
%       SHIFTS : for multiresolution (eye(length(DIMS_K)))
%       DEBUG : bool to prompt procedure state (false)
%       STABILDIAG : bool to plot the stabilization diagram (false)
%       MAC : correlation coefficient to link the poles by MAC values (0)
%           - MAC == 0 : no MAC value computed
%           - 0 < MAC < 1 : poles with corr > MAC are linked
%
%   OUT structure : 
%    % Sizes
%       Lp : dimensions of the signal over the isophase surface (size : [1 length(DIMS_P)])
%       Lk : dimensions of the signal over the wavevectors (size : [1 length(DIMS_K)])
%       Kk : number of rows of the signal matrix for each signal dimension (size : [1 length(DIMS_K)])
%       Mk : number of columns of the signal matrix for each signal dimension (size : [1 length(DIMS_K)])
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
        K0 ;
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
    
    % Reshaping the signal
        Signal = permute(Signal,[DIMS_P DIMS_K]) ; % Sort dimensions
        Signal = reshape(Signal,[prod(Lp) prod(Lk)]) ;
        [~,permDims] = sort([DIMS_P DIMS_K]) ; % Used to rechape the signal model
        
    % Output Initialization
        K = zeros([max(R0) length(DIMS_K)]) ;
        dK = zeros([max(R0) length(DIMS_K)]) ;
        U = zeros([max(R0) length(DIMS_P)+1]) ;
        SignalModel = zeros(size(Signal)) ;
        ESTER = zeros([size(SHIFTS,1) length(R0)]) ;
        
        
        
        
% HANKEL MATRIX BUILDING
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Hankel Matrix Initialization : ') ; lastTime = tic ; end
    
    % Sub-Hankel matrices shapes (Kk rows, Mk columns)
        DECIM_K = DECIM(DIMS_K) ;
        Kk = zeros([1 length(DIMS_K)]) ; % max size of (signal+noise) subspace
        Mk = zeros([1 length(DIMS_K)]) ;
        isCOS = false([1 length(DIMS_K)]) ; % IS COSINUS SEARCHED ?
        for d = 1:length(DIMS_K)
            switch lower(FUNC(d,:))
                case 'exp'
                    if isempty(K0) 
                        Kk(d) = floor((Lk(d)+DECIM_K(d))./(1+DECIM_K(d))) ;
                    else
                        Kk(d) = K0(d) ;
                    end
                    Mk(d) = Lk(d) - DECIM_K(d).*(Kk(d)-1) ;
                case 'cos'
                    if isempty(K0) 
                        Kk(d) = floor((Lk(d)+DECIM_K(d))/(2*DECIM_K(d)+1)) ;
                    else
                        Kk(d) = K0(d) ;
                    end
                    Mk(d) = Lk(d) - DECIM_K(d)*(2*Kk(d)-1) ;
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
                case 'cos'
                    subIndH{d,1} = ((Kk(d)-1)*DECIM_K(d)+1:-DECIM_K(d):1)'*ones(1,Mk(d)) + ones(Kk(d),1)*(0:Mk(d)-1) ;
                    subIndH{d,2} = (Kk(d)*DECIM_K(d):DECIM_K(d):DECIM_K(d)*(2*Kk(d)-1))'*ones(1,Mk(d)) + ones(Kk(d),1)*(1:Mk(d)) ;
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
            
        % Complete indices
            indHbH = [] ;
            n_indHbH = [] ;
            indP = [] ;
            indicesHbH ;
            
        
        
        
% AUTOCOVARIANCE MATRIX
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Computing of Css : ') ; lastTime = tic ; end
  
        Css = buildCss ;
        
        
        
% SIGNAL SUBSPACE ESTIMATION        
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Eigendecomposition of Css : ') ; lastTime = tic ; end
    % ESTIMATION
        if isempty(W0) % re compute completely the signal subspace
            switch SOLVER
                case 'eig'
                    [W,lambda] = eig(Css,'vector') ;
                case 'eigs'
                    [W,lambda] = eigs(Css,max(R0),'lm') ;
                    lambda = diag(lambda) ;
                case 'svds'
                    H = buildHss(Signal) ;
                    [Us,Sigma,Vs] = svds(H,max(R0)) ;
                    lambda = diag(Sigma).^2 ;
                    W = Us ;
                case 'mysvd'
                    H = buildHss(Signal) ;
                    if  size(H,1) <= size(H,2)
                        C = H*H';
                        [Us,Lambda] = eigs(C,max(R0));
                        [lambda,ix] = sort(abs(diag(Lambda)),'descend');
                        Us = Us(:,ix);   
                        Vs = H'*Us;
                        sigma = sqrt(lambda);
                        Vs = bsxfun(@(x,c)x./c, Vs, sigma');
                        Sigma = diag(sigma);
                    else
                        C = H'*H; 
                        [Vs,Lambda] = eigs(C,max(R0));
                        [lambda,ix] = sort(abs(diag(Lambda)),'descend');
                        Vs = Vs(:,ix);    
                        Us = H*Vs;
                        sigma = sqrt(lambda);
                        Us = bsxfun(@(x,c)x./c, Us, sigma');
                        Sigma = diag(sigma);
                    end
                    W = Us ;
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
                    CRIT = [MDL(R0) ; 0] ;
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
    extractPoles(R) ;
        
    
% ESTIMATE AMPLITUDES
    if COMPUTE_U || COMPUTE_dK
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Amplitude Estimation : ') ; lastTime = tic ; end
        computeU ;
        A ;
    end
    
    
% ESTIMATE UNCERTAINTIES
    if COMPUTE_dK 
        if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Uncertainty Estimation : ') ; lastTime = tic ; end
        computeUncertainties ;
    end
        
    
% OUTPUT
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('   Output : ') ; lastTime = tic ; end
    % Results
        OUT.K = K ;
        OUT.dK = dK ;
        OUT.U = U ;
        if isempty(DIMS_P)
            OUT.SignalModel = reshape(SignalModel,Lk) ;
        else
            OUT.SignalModel = permute(reshape(SignalModel,[Lp Lk]),permDims) ;
        end
    % Sizes
        OUT.Lp = Lp ;
        OUT.Lk = Lk ;
        OUT.Kk = Kk ;
        OUT.Mk = Mk ;
    % Subspace
        OUT.Css = Css ;
        OUT.lambda = lambda ;
        OUT.W = W ;
    

% END OF THE PROCEDURE
    if(DEBUG) ; display(['       ',num2str(toc(lastTime),3), ' sec']) ; display('---------------------------------') ; end


    
    
    
    
    
% ===================================================================================================================    
% FUNCTIONS FOR ESPRIT
% ===================================================================================================================

% HANKEL-BLOCK-HANKEL MATRIX INDICES
    function indicesHbH
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
    end

% BUILD THE CSS MATRIX
    function Css = buildCss
        % Build
            if any(Kk~=Lk) % Mk~=1, spatial smoothing
                % Sum the elementary Css matrices by points of isosurface
                    Css = zeros(prod(Kk),prod(Kk)) ;
                    for p = 1:length(indP)
                        sig = Signal(indP(p),:) ;
                        SIG = sig(indHbH{1}) ;
                        for i = 2:n_indHbH
                            SIG = SIG + sig(indHbH{i}) ;
                        end
                        SIG = SIG/n_indHbH ;
                        Css = Css + SIG*SIG' ;
                    end
            else % Mk==1, no spatial smoothing (multiple snapshots case)
                SIG = Signal.' ;
                Css = SIG*SIG' ;
            end
        % Normalize
            Css = Css/(prod(Mk)*length(indP)) ;
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
                Hss = Data.' ;
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
                    I = eye(prod(Kk)) ;
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
            shiftsCOS = logical(repmat(isCOS(:)',[size(SHIFTS,1) 1])) ;
            %shiftsCOS = logical(repmat(isCOS(:),[1 1])) ;
            if 0
                K = zeros(size(Z)) ;
                K(~shiftsCOS,:) = log(Z(~shiftsCOS,:))/1i ; % FUNC = 'EXP' ;
                K(shiftsCOS,:) = acos(Z(shiftsCOS,:)) ; % FUNC = 'COS' ;
            else
                K = log(Z)/1i ;
            end
        % Wavevectors in the cartesian basis
            K = (SHIFTS*diag(DECIM_K))\(K) ;
    end



% VANDERMONDE MATRIX
    function V = buildVandermonde(L) % Here, decimation factors are no longer taken into account
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
            V = buildVandermonde(Lk) ;
        % Shift for the conditionning
            v0 = V(1,:) ;
            V = V*diag(1./v0) ;
        % Signal reshaping
            S = Signal.' ;
        % Amplitude estimation
            A = V\S ;
        % Vandermonde Shift conpensation
            A = diag(v0)*A ;
        % Reshaping
            U = reshape(A.',[Lp , size(K,2)]) ;
        % "Pure" signal reconstruction
            SignalModel = (V*A).' ;
    end


% UNCERTAINTY ESTIMATION dK
%     function computeUncertainties
%         % Init
%             dK = zeros(size(K)) ;
%         % Delta Signal
%             dS = Signal-SignalModel ;
%         % Delta signal matrix
%             dH = buildHss(dS) ;
%         % Pure signal matrix
%             H = buildHss(SignalModel) ;
%         % Signal subspace of the noiseless signal
%             [Us,omega] = eig(H*H'/prod(Mk)/length(indP),'vector');
%             [~,ind] = sort(omega,'descend') ;
%             ind = ind(1:R0(R)) ;
%             Us = Us(:,ind) ;
%             omega = omega(ind) ;
%         % Signal subspace perturbation
%             dUs = (eye(size(dH,1))-Us*Us')*dH*pinv(H)*Us ;%
%             %dUs = W-Us ;%
%         % Uncertainties
%             T = 1 ;
%             invT = 1 ;
% %             [T,~] = extractPoles(R) ;
% %             invT = inv(T) ;
%             shifts = eye(length(DIMS_K)) ; % /!\ SHIFTS IS DEFAULT HERE !
%             for n = 1:size(K,1)
%                 [~,~,~,Jup,Jdwn] = selectMatrices(shifts(n,:)) ;
%                 Fn = (Jup*Us)\(Jdwn*Us) ;
%                 dFn = Fn*(Jdwn*Us)\(Jdwn*dUs) - (Jup*Us)\(Jup*dUs)*Fn ;
%                 %dFn = (Jup*W)\(Jdwn*W) - (Jup*Us)\(Jdwn*Us) ;
%                 for o = 1:size(K,2)
%                     darn = invT(:,o).'*dFn*T(:,o) ;
%                     dK(n,o) = abs(darn/exp(1i*K(n,o)))^2 ;
%                     %dK(s,r) = abs(darn/K(s,r))^2 ;
%                 end
%             end 
%     end


% UNCERTAINTY ESTIMATION dK
    function computeUncertainties
        % Init
            dK = zeros(size(K)) ;
        % Delta Signal
            dS = Signal-SignalModel ;
        % Delta signal matrix
            dH = buildHss(dS)*sqrt(R0(R))/prod(DECIM_K) ; % %
        % Partial Vandermonde matrices
            V = buildVandermonde(Lk) ;
            P = V(indHbH{1}(:,1),:) ;
            Q = V(indHbH{1}(1,:),:) ;
            for i = 2:n_indHbH % If cosinuses
                P = P + V(indHbH{2}(:,1),:) ; 
            end
            %P = P*(norm(P)) ;
            %Q = Q*(norm(Q)) ;
        % Complete right-Vandermonde Matrix
            QQ = zeros(prod(Mk)*length(indP),R0(R)) ;
            for p = 1:length(indP)
                QQ((1:prod(Mk))+prod(Mk)*(p-1),:) = Q*diag(A(:,indP(p))) ;
            end
        % Uncertainties
            shifts = eye(length(DIMS_K)) ; % /!\ SHIFTS IS DEFAULT HERE !
            for s = 1:size(K,1)
                [~,~,~,Jup,Jdwn] = selectMatrices(shifts(s,:)) ;
                for r = 1:size(K,2)
                    br = [zeros(r-1,1) ; 1 ; zeros(size(K,2)-r,1)] ;
                    arn = exp(1i*K(s,r)) ;
                    vrn = br'*((Jup*P)\(Jdwn-arn*Jup)) ;
                    xr = (QQ.')\br ;
                    dK(s,r) = abs(vrn*dH*xr/arn)^2/(prod(Mk)*length(indP)*max(1,prod(Mk./Kk))) ;
                end
            end 
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
                        case 'K0'
                            K0 = Value ;
                            paramSet(15) = true ;
                        case 'W0'
                            W0 = Value ;
                            paramSet(16) = true ;
                        case 'COMPUTE_DK'
                            COMPUTE_dK = Value ;
                            paramSet(17) = true ;
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
            if ~paramSet(14) ; COMPUTE_U = true ; end
            if ~paramSet(15) ; K0 = [] ; end
            if ~paramSet(16) ; W0 = [] ; end
            if ~paramSet(17) ; COMPUTE_dK = true ; end
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
                    plot3(abs(real(K_MAClinks')),OR_MAClinks',abs(imag(K_MAClinks')./real(K_MAClinks')),'-k','linewidth',.5)
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
                % Limits
                    %set([stab.axMean stab.axPoles],'xlim',[1/prod(Lk)*2*pi/2 pi]) ;
                    w = abs(real(Kstab(:))) ;
                    set([stab.axMean stab.axPoles],'xlim',[min(w(w~=0)) max(w)].*[0.9 1.1]) ;
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