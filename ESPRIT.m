classdef ESPRIT < handle
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

%% INPUT (SET ACCESSIBLE) PARAMETERS
properties
    Signal
    DIMS_K(1,:) 
    M_L(1,:) = 1/2
    DECIM(1,:)
    R0(1,:) 
    W0
    SOLVER char = 'tensor'
    CRITERION char = 'MDL'
    CRIT_THRS(1,1) = 1
    CHOICE char = 'auto'
    FUNC cell
    SHIFTS 
    FIT(1,:) char = 'TLS'
    DEBUG(1,1) logical = false
    COMPUTE_U(1,1) logical = false
    STABILDIAG(1,1) logical = false
    SIGNAL_MODEL(1,1) logical = false
    COMPUTE_dU(1,1) logical = false
    COMPUTE_dK(1,1) logical = false
    MAC(1,1) logical = false
    UserData = [] % for user-custom data
end

%% OUTPUT PROPERTIES (FOR COMPATIBILITY REASONS)
properties
    K,U,PHI,dK,dU
end

%% CONSTRUCTOR/DESTRUCTOR
methods
    function this = ESPRIT(varargin)
    % Class Constructor
        if nargin==0 ; return ; end
    % Depending on the number of numeric pre-inputs...
        nPreNumeric = find(cellfun(@ischar,[varargin 'end']),1,'first')-1 ; 
    % ESPRIT(Signal,...)
        if nPreNumeric>=1 ; this.Signal = varargin{1} ; end 
    % ESPRIT(Signal,DIMS_K,...)
        if nPreNumeric>=2 ; this.DIMS_K = varargin{2} ; end
    % ESPRIT(Signal,DIMS_K,R0,...)
        if nPreNumeric>=3 ; this.R0 = varargin{3} ; end
    % Set properties
        varargin = varargin(nPreNumeric+1:end) ;
        this.setProperties(varargin{:}) ;
    % If no signal is available, do not go further
        if isempty(this.Signal) ; return ; end
    % Extract the wavevectors
        this.K = this.extractPoles ; 
    % Extract the modes
        if this.COMPUTE_U
            V = this.vandermondeMatrix(this.K) ;
            A = this.extractAmplitudes(V) ;
            [this.U,this.PHI] = this.extractModes(A) ;
        end
    % Compute uncertainties
        if this.COMPUTE_dK || this.COMPUTE_dU
            [this.dK,this.dU] = this.computeUncertainties ;
        end
    end
    
    function delete(this)
    % Class Destructor
    end
end

%% SET PROPERTIES
methods
    function setProperties(this,varargin)
    % Set the properties with a given argument array
        if mod(numel(varargin),2) ; error('Wrong number of arguments provided') ; end
        propNames = varargin(1:2:end-1) ;
        propValues = varargin(2:2:end) ;
        for pp=1:numel(propNames)
            this.(propNames{pp}) = propValues{pp} ;
        end
    % Set the other properties
        this.setDefaultProperties ;
    end
    
    function setDefaultProperties(this)
    % Set the default value of dependent properties
        if isempty(this.Signal) ; return ; end
        if isempty(this.DIMS_K) ; this.DIMS_K = find(size(this.Signal)>1,1,'last') ; end 
        if isempty(this.DECIM) ; this.DECIM = ones(1,ndims(this.Signal)) ; end  
        if isempty(this.R0) ; this.R0 = 1:floor(min(arrayfun(@(d)size(this.Signal,d),this.DIMS_K)/2)) ; end
        if isempty(this.FUNC) ; this.FUNC = repmat({'exp'},[numel(this.DIMS_K) 1]) ; end 
        if isempty(this.SHIFTS) ; this.SHIFTS = eye(length(this.DIMS_K)) ; end   
    end
end


%% SIZE INFORMATION
properties (SetAccess=protected,Dependent)
    DIMS_P(1,:) % dimensions of the signal corresponding to amplitudes/modes
    Lk(1,:) % total length of the signal in the wavevector dimensions
    Lp(1,:) % total length of the signal in the amplitudes/modes dimensions
    DECIM_K(1,:) % decimation factors on wavevector dimensions
    DECIM_P(1,:) % decimation factors on amplitudes/modesr dimensions
end
methods
    function out = get.DIMS_P(this) 
        out = setdiff(1:ndims(this.Signal),this.DIMS_K) ; 
        if isempty(out) ; out = ndims(this.Signal) + 1 ; end
    end
    function out = get.Lk(this) ; out = size(this.Signal,this.DIMS_K) ; end
    function out = get.Lp(this) ; out = size(this.Signal,this.DIMS_P) ; end
    function out = get.DECIM_K(this) ; out = this.DECIM(this.DIMS_K) ; end
    function out = get.DECIM_P(this) 
        dec = [this.DECIM 1] ;
        out = dec(this.DIMS_P) ; 
    end
% Which dimension corresponds to an exponential/cosine component ? 
    function out = isCOS(this) ; out = ismember(this.FUNC(:)','cos') ; end
    function out = isEXP(this) ; out = ismember(this.FUNC(:)','exp') ; end
end


%% SPATIAL SMOOTHING LENGTHS
properties (SetAccess=protected,Dependent)
    Mk(1,:) % Spatial smoothing snapshots
    Kk(1,:) % Snapshot size
end
methods
    function m_l = spatialSmoothing(this)
    % Correct spatial smoothing criterion
        m_l = this.M_L ;
        m_l = max(m_l,1./this.Lk) ; % avoids Mk=0
        m_l = min(m_l,1-1./this.Lk) ; % avoids Kk=0
    end
    
    function k = get.Kk(this) 
    % Snapshot size: max size of (signal+noise) subspace
        k = floor((1-this.spatialSmoothing).*this.Lk./this.DECIM_K + 1) ;        
    % When Decimation, force the spatial smoothing to use all data
        decimAndNoSS = this.Lk-(k-1).*this.DECIM_K<this.DECIM_K ;
        k(decimAndNoSS) = k(decimAndNoSS)-1 ;
    % Cosines: use half the data
        k(this.isCOS) = floor(k(this.isCOS)/2) ;
    end
    
    function out = get.Mk(this) 
    % Number of spatial smoothing snapshots
        k = this.Kk ;
        k(this.isCOS) = k(this.isCOS)*2 ;
        out = this.Lk - this.DECIM_K.*(k-1) ;
    end
end


%% HANKEL MATRIX INDICES
methods
    function nh = nSubHankelMat(this,fun)
        if nargin<2 ; fun = this.FUNC ; end
        nh = any(ismember(fun,'cos'))+1 ; 
    end
    
    function IND = hankelSubIndices(this,k,m,s,fun)
    % Indices of hankel matrices associated to individual signal dimensions
        if nargin<2 ; k = this.Kk ; end
        if nargin<3 ; m = this.Mk ; end
        if nargin<4 ; s = this.DECIM_K ; end
        if nargin<5 ; fun = this.FUNC ; end
    % Sub indices
        nD = numel(k) ;
        IND = cell(nD,this.nSubHankelMat(fun)) ; % cell array of [Kk(d) Mk(d)] indices
        for d = 1:nD
            switch lower(fun{d})
                case 'exp'
                    IND{d,1} = s(d)*(0:k(d)-1)' + (0:m(d)-1) + 1 ;
                case 'cos'
                    IND{d,1} = s(d)*(k(d)-1:-1:0)' + (0:m(d)-1) + 1 ;
                    IND{d,2} = s(d)*(k(d):1:2*k(d)-1)' + (1:m(d)) + 1 ;
            end
        end
    end
    
    function IND = blockHankelSubIndices(this,k,m,s,fun)
    % Block Hankel matrix sub-indices (with kronecker struture) associated to each signal dimension
        if nargin<2 ; k = this.Kk ; end
        if nargin<3 ; m = this.Mk ; end
        if nargin<4 ; s = this.DECIM_K ; end
        if nargin<5 ; fun = this.FUNC ; end
    % Sub indices
        HI = this.hankelSubIndices(k,m,s,fun) ;
        ck = cumprod([1 k(1:end-1)]) ; ckr = cumprod([k(2:end) 1],'reverse') ;
        cm = cumprod([1 m(1:end-1)]) ; cmr = cumprod([m(2:end) 1],'reverse') ;
    % All Indices
        nD = numel(k) ;
        IND = cell(nD,size(HI,2)) ; % cell array of [prod(Kk) prod(Mk)] indices
        for d = 1:nD
            switch lower(fun{d})
                case 'exp'
                    IND{d,1} = repelem(HI{d,1},ck(d),cm(d)) ;
                    IND{d,1} = repmat(IND{d,1},ckr(d),cmr(d)) ;
                case 'cos'
                    for i = 1:2
                        IND{d,i} = repelem(HI{d,i},ck(d),cm(d)) ;
                        IND{d,i} = repmat(IND{d,i},ckr(d),cmr(d)) ;
                    end
            end
        end
    end
    
    function IND = blockHankelIndices(this,k,m,s,fun)
    % Block Hankel matrix LINEAR indices
        if nargin<2 ; k = this.Kk ; end
        if nargin<3 ; m = this.Mk ; end
        if nargin<4 ; s = this.DECIM_K ; end
        if nargin<5 ; fun = this.FUNC ; end
    % Cummulative numer of elements before the dimension
        lk = m + (k-1).*s ;
        clk = cumprod(lk) ;
    % First dimension
        bHSI = this.blockHankelSubIndices(k,m,s,fun) ;
        IND = bHSI(1,:) ;
    % Other dimensions
        nD = numel(k) ;
        if nD>1
        % EXPONENTIALS
            for d = 2:nD
                IND{1} = IND{1} + (bHSI{d}-1)*clk(d-1) ;
            end
        % COSINES (IF NEEDED)
            ic = ismember(fun,'cos') ;
            if any(ic)
                for d = 2:nD
                    IND{2} = IND{2} + (bHSI{d,1+ic(d)}-1)*clk(d-1) ;
                end
            end
        end
    end
end


%% INDICES FOR SHIFT INVARIANCE
methods
    function [indUp,indDwn1,indDwn2] = shiftIndices(this,shift)
    % Build the vector of selection indices corresponding to a shift invariance
        shift_dims = find(shift) ;
        isShiftCOS = any(this.isCOS.*shift) ; % Cosinus searched in the shift ?
    % Subindices of the first colon of the H-b-H matrix
        k = this.Kk ;
        dec = this.DECIM_K ;
        hsi = this.blockHankelSubIndices(k,k*0+1,dec) ;
    % Indices Up
        indUp = true(size(hsi{1})) ; 
        if isShiftCOS % Cosinus searched in the shift
            indTrUp = true(size(hsi{1})) ; 
        end
        for d = shift_dims
            switch lower(this.FUNC{d})
                case 'exp'
                    if shift(d)>0 ; indUp = indUp & (hsi{d,1}<=(k(d)-shift(d))*dec(d)) ; end
                    if shift(d)<0 ; indUp = indUp & (hsi{d,1}>-shift(d)*dec(d)) ; end
                case 'cos'
                    indTrUp = indTrUp & (hsi{d,1}<=(k(d)-shift(d))*dec(d)) ;
                    indTrUp = indTrUp & (hsi{d,1}>shift(d)*dec(d)) ;
            end
        end
    % Indices Down
        indDwn = true(size(hsi{1})) ; 
        if isShiftCOS % Cosinus searched in the shift
            indTrDwn1 = true(size(hsi{1})) ; 
            indTrDwn2 = true(size(hsi{1})) ; 
        end
        for d = shift_dims
            switch lower(this.FUNC{d})
                case 'exp'
                    if shift(d)>0 ; indDwn = indDwn & (hsi{d,1}>shift(d)*dec(d)) ; end
                    if shift(d)<0 ; indDwn = indDwn & (hsi{d,1}<=(k(d)+shift(d))*dec(d)) ; end
                case 'cos'
                    indTrDwn1 = indTrDwn1 & (hsi{d,1}<=(k(d)-2*abs(shift(d)))*dec(d)) ;
                    indTrDwn2 = indTrDwn2 & (hsi{d,1}>2*abs(shift(d))*dec(d)) ;
            end
        end
    % Combine Indices
        if isShiftCOS % Cosinus searched in the shift
            indUp = indUp & indTrUp ;
            indDwn1 = indDwn & indTrDwn1 ;
            indDwn2 = indDwn & indTrDwn2 ;
        else % Only exponentials in the shift
            indDwn1 = indDwn ;
            indDwn2 = [] ;
        end
    end
    
    function [Jup,Jdwn] = selectionMatrices(this,shift)
    % Selection matrices corresponding to a shift
        I = speye(prod(this.Kk)) ;
        [indUp,indDwn1,indDwn2] = shiftIndices(this,shift) ;
        % Jup
            Jup = I(indUp,:) ;
        % Jdwn
            Jdwn = I(indDwn1,:) ;
            if any(this.isCOS.*shift) % Cosinus searched in the shift
                Jdwn = ( Jdwn + I(indDwn2,:) )/2 ;
            end
    end
end


%% SIGNAL RESHAPING
methods
    function S = reshapedSignal(this,S)
    % Return the reshaped version of the signal
        if nargin<2 ; S = this.Signal ; end
        S = permute(S,[this.DIMS_K this.DIMS_P]) ;
        S = reshape(S,[prod(this.Lk) prod(this.Lp)]) ;
    end
    
    function d = permDims(this)
    % Permutation dimensions of the Signal<->reshapedSignal
        d = sort([this.DIMS_K this.DIMS_P]) ;
    end
end


%% COVARIANCE
methods
    function indP = modeIndices(this)
    % Return the points included in the hankel signal matrix
    % Decimation of the isosurface/modes/amplitudes points
        dP = this.DECIM_P ;
        if any(dP~=1) % SOME DECIMATION IS NEEDED
            lp = this.Lp ;
            clp = cumprod(lp) ;
            indP = 1:dP(1):lp(1) ;
            for d = 2:length(dP)
                indP = indP(:) + clp(d-1)*(0:dP(d):lp(d)-1) ;
            end
            indP = reshape(indP,1,[]) ;
        else % NO DECIMATION
            indP = 1:prod(size(this.Signal,this.DIMS_P)) ;
        end
    end
    
    function Hs = signalMatrix(this,S,iHbH)
    % Return the H-b-H signal matrix 
        if nargin<2 ; S = this.Signal ; end
        if nargin<3 ; iHbH = this.blockHankelIndices ; end
    % Take the signal decimated on points
        S = this.reshapedSignal(S) ;
        S = S(:,this.modeIndices) ;
    % Build the Hankel matrix
        if any(this.Mk~=1) % Mk~=1, spatial smoothing
            Hs = S(iHbH{1},:) ;
            for i = 2:numel(iHbH) ; Hs = Hs + S(iHbH{i},:) ; end
        else % Mk==1, no spatial smoothing (multiple snapshots case)
            Hs = S ;
        end
    % Reshape size: [prod(Kk) prod(Mk)*numel(indP)]
        Hs = reshape(Hs,prod(this.Kk),[]) ;
    end
    
    function Css = covarianceMatrix(this,S,iHbH,method)
    % Compute the covariance matrix of the signal
        if nargin<2 ; S = this.Signal ; end
        if nargin<3 ; iHbH = this.blockHankelIndices ; end
        if nargin<4 ; method = 'prod' ; end
        k = this.Kk ; m = this.Mk ;
    % Take the signal decimated on points
        S = this.reshapedSignal(S) ;
        S = S(:,this.modeIndices) ;
    % Compute Css
        if any(m~=1) % Mk~=1, spatial smoothing
        % H-b-H indices
            S = S/numel(iHbH) ;
        % Covariance Css=Hs*Hs'
            switch method
                case 'prod' % Sum the elementary Css matrices by points of isosurface
                    Css = zeros(prod(k)) ;
                    %t = tic ; 
                    for p = 1:size(S,2)
                        H = S(iHbH{1},p) ;
                        for i = 2:numel(iHbH)
                            H = H + S(iHbH{i},p) ;
                        end
                        H = reshape(H,prod(k),[]) ;
                        Css = Css + H*H' ;
                        %if length(indP)>1 && p==1 && this.DEBUG ; disp(['       estim: ',num2str(toc(t)*length(indP)),' secs']) ; end
                    end
                case 'conv' % use convolution
                    error('Convolution is not implemented yet.') ;
            end
        else % Mk==1, no spatial smoothing (multiple snapshots case)
            Css = S*S' ;
        end
    end
end


%% SIGNAL ORDER CRITERION
methods
    function R = signalOrder(this,R0,W,lambda,criterion,thrs)
    % Return the signal order given by a criterion
        if nargin<2 ; R0 = this.R0 ; end
        if numel(R0)==1 ; R = R0 ; return ; end
        if nargin<4 ; [W,lambda] = this.signalSubspace(max(R0)+1) ; end
        if nargin<5 ; criterion = this.CRITERION ; end
        if nargin<6 ; thrs = this.CRIT_THRS ; end
    % Compute the criterion
        switch upper(criterion)
            case 'MDL'
                crit = this.mdl(lambda) ;
            case 'ESTER'
                crit = this.ester(W,this.SHIFTS) ;
            case 'SAMOS'
                crit = this.samos(W,this.SHIFTS) ;
        end
    % Mean over shifts
        crit = prod(crit,1) ;
    % Select the signal order
        crit = crit(R0) ;
        r = find(crit>=max(crit)*thrs,1,'last') ;
        R = R0(r) ;
    end
    
    function m = mdl(this,lambda)
    % Minimum Definition Length criterion
        if nargin<2 ; [~,lambda] = signalSubspace(this) ; end
        %m = -diff(lambda)./lambda(1:end-1) ;
        m = lambda(1:end-1)./lambda(2:end) ;
        m = m(:)' ;
    end
    
    function err = ester(this,W,shifts)
    % ESTimation of ERror
        if nargin<2 ; W = this.signalSubspace ; end
        if nargin<3 ; shifts = this.SHIFTS ; end
        err = NaN(size(shifts,1),size(W,2)) ;
        for r = 1:size(W,2)
            Wr = W(:,1:r) ;
            for s = 1:size(shifts,1)
                [Wup,Wdwn,F] = this.spectralMatrix(shifts(s,:),Wr) ;
                err(s,r) = 1/norm(Wup*F-Wdwn) ;
            end
        end
    end

    function err = samos(this,W,shifts)
    % Subspace Automatic Model Order Selection
        if nargin<2 ; W = this.signalSubspace ; end
        if nargin<3 ; shifts = this.SHIFTS ; end
        err = NaN(size(shifts,1),size(W,2)) ;
        for r = 1:size(W,2)
            Wr = W(:,1:r) ;
            for s = 1:size(shifts,1)
                [Wup,Wdwn] = this.spectralMatrix(shifts(s,:),Wr) ;
                S = [Wup Wdwn] ;
                g = sqrt(real(sort(eig(S'*S,'vector')))) ;
                err(s,r) = r/sum(g(1:r)) ;
            end
        end
    end
end


%% SIGNAL SUBSPACE
methods
    function [W,lambda] = signalSubspace(this,R,method)
    % Estimate the signal subspace
        if nargin<2 ; R = max(this.R0)+1 ; end
        if nargin<3 ; method = this.SOLVER ; end
    % Pre-computations
        if ~ismember(method,{'tensor'})
            Css = this.covarianceMatrix ;
        end
    % Adjust signal order candidates
        k = this.Kk ; m = this.Mk ; indP = this.modeIndices ;
        Rmax = min(prod(k),prod(m)*numel(indP)) ;
        Rmax = Rmax - 1- max(sum(abs(this.SHIFTS),2)) ;
        R = min(R,Rmax) ;
    % ESTIMATION
        switch method
            case 'eig' % Re-compute the COMPLETE decomposition
                [W,lambda] = eig(Css,'vector') ;
            case 'eigs'  % Re-compute the signal subspace ONLY
                [W,lambda] = eigs(Css,R,'lm') ;
                lambda = diag(lambda) ;
            case 'tensor' % High-order tensor decomposition
            % Build the signal tensor
                Hs = this.signalMatrix ;
                nD = numel(k) ;
                A = reshape(Hs,[k prod(m)*numel(indP)]) ;
            % High-Order SVD in the nD dimensions
                [UU,S,ll] = this.hosvd(A,1:nD+1,R) ;
            % Build the signal space
                W = this.tuckerprod(1:nD,S,UU(1:nD)) ;
                W = reshape(W,[prod(k) R]) ;
                lambda = ll{nD+1} ;
            case 'follow' % Approximate with one QR iteration (Badeau-style)
                Cxy = Css*this.W0 ;
                [W,~] = qr(Cxy,0) ;
                lambda = diag(W'*Css*W) ;
        end
    % SORT EIGENVALUES IN DESCENDING ORDER
        [~,ind] = sort(lambda,'descend') ;
        lambda = lambda(ind) ;
        W = W(:,ind) ;
    % KEEP ONLY THE HIGHEST E.V
        W = W(:,1:R) ; 
        lambda = lambda(1:R) ;
    end
end
methods (Static)
    function A = unfold(A,d)
    % Return the mode-d unfolding of a tensor A
        A = permute(A,[d 1:d-1 d+1:ndims(A)]) ;
        A = reshape(A,size(A,1),[]) ;
    end
    
    function A = refold(A,d,sz)
    % Return the mode-d folded tensor
        foldDims = [d 1:d-1 d+1:numel(sz)] ;
        A = reshape(A,sz(foldDims)) ;
        A = permute(A,[2:d 1 d+1:numel(sz)]) ;
    end
    
    function B = tuckerprod(d,A,U)
    % Return the d-mode Tucker product of a tensor with a matrix 
    % Recursive multi-product ? (U has to be a cell array of matrices)
        if numel(d)>1
            A = ESPRIT.tuckerprod(d(1:end-1),A,U(1:end-1)) ; 
        end
    % Select only the last dimension
        d = d(end) ;
        if iscell(U) ; U = U{end} ; end
    % Size information
        szB = size(A) ; 
        szB(d) = size(U,1) ;
    % Tucker product
        B = U*ESPRIT.unfold(A,d) ;
        B = ESPRIT.refold(B,d,szB) ;
    end
    
    function [U,S,lambda] = hosvd(A,dims,R)
    % High-order (truncated) SVD decomposition of a tensor
    % Dimensions along which to apply the SVD
        sz = size(A) ; nD = ndims(A) ;
        if nargin<2 || isempty(dims) ; dims = 1:nD ; end
    % R is the truncation order
        if nargin<3 || isempty(R) ; R = size(A,dims) ; end
        R = R(:)'.*ones(1,numel(dims)) ;
    % Perform the decomposition (interlaced implementation)
        U = cell(1,nD) ; lambda = cell(1,nD) ;
        S = A ;
        for dd = 1:numel(dims)
            d = dims(dd) ;
        % mode-d unfolding of the tensor
            Ad = ESPRIT.unfold(S,d) ;
        % covariance
            leftEV = size(Ad,1)<=size(Ad,2) ;
            if leftEV % Usual covariance, compute left eigenvectors U
                Cdd = Ad*Ad' ;
            else % Compute the right-eigenvectors V
                Cdd = Ad'*Ad ; 
            end
        % eigen space associated to the mode-d unfold
            %[U{d},lambda{d}] = eigs(Cdd,min(R(dd),sz(d)),'lm') ;
            [U{d},lambda{d}] = eig(Cdd,'vector') ;
            lambda{d} = sqrt(lambda{d}) ; % eigenvalues of Ad²
        % left-eigenvectors
            if ~leftEV % Ad = U*diag(l)*V' -> U = Ad*V*diag(1./l) ;
                U{d} = Ad*U{d}*diag(1./lambda{d}) ;
            end
        % Select the R dominant components if needed
            if R(dd)<sz(d) % truncate ?
                [~,is] = sort(lambda{d},'descend') ;
                is = is(1:R(dd)) ;
                lambda{d} = lambda{d}(is) ;
                U{d} = U{d}(:,is) ;
                sz(d) = R(dd) ;
            end
        % reduce the signal tensor (at the end, A is the core tensor)
            S = ESPRIT.refold(U{d}'*Ad,d,sz) ;
        end
    end
end


%% SHIFT INVARIANCE
methods
    function [Wup,Wdwn,F] = spectralMatrix(this,shift,W)
    % Return the spectral matrix associated to a given shift
        if nargin<3 ; W = this.signalSubspace ; end
        R = size(W,2) ;
    % Get shift indices
        [indUp,indDwn1,indDwn2] = this.shiftIndices(shift) ;
    % W_up and W_down construction
        if any(this.isCOS.*shift) 
            Wup = W(indUp,:) ;
            Wdwn = ( W(indDwn1,:) + W(indDwn2,:) )/2 ;
        else
            Wup = W(indUp,:) ;
            Wdwn = W(indDwn1,:) ;
        end
    % Only the Wup and Wdwn matrices needed ?
        if nargout==2 ; return ; end
    % Shift invariance evaluation
        switch this.FIT
            case 'LS' % Least-Squares
                F = Wup\Wdwn ;
            case 'TLS' % Total Least-Squares
                [E,~] = svd([Wup' ; Wdwn']*[Wup Wdwn],0) ;
                F = -E(1:R,R+(1:R))/E(R+(1:R),R+(1:R)) ;
            case 'TLS2' % TLS with less correlation between samples (experimental)
                L_2 = floor(this.Lk./2) ;
                Wup = W(1:2:L_2*2,:) ;
                Wdwn = W(2:2:L_2*2,:) ;
                [E,~] = svd([Wup' ; Wdwn']*[Wup Wdwn],0) ;
                F = -E(1:R,R+(1:R))/E(R+(1:R),R+(1:R)) ;
        end
    end
end


%% POLE EXTRACTION
methods
    function [K,T,Beta] = extractPoles(this,W,R,shifts)
    % Extract the wavevectors from the signal
        if nargin<2 ; [W,lambda] = this.signalSubspace ; end
        if nargin<3 ; R = this.signalOrder(this.R0,W,lambda) ; end
        if nargin<4 ; shifts = this.SHIFTS ; end
        nS = size(shifts,1) ;
    % Evaluation of all shift invariances
        W = W(:,1:R) ;
        F = zeros(R,R,nS) ;
        for s = 1:nS
            [~,~,F(:,:,s)] = this.spectralMatrix(shifts(s,:),W) ;
        end
    % Evaluation of PHI
        if nS==1 % One shift, simple diagonalization
            [T,PI] = eig(F,'vector') ;
            Beta = 1 ;
        else % More than one shift, joint diagonalization 
        % with a linear combination of F's
        % Coefficients of the linear combination
            Beta = 1 + 1.2.^(0:nS-1) ;
            Beta = reshape(Beta/sum(Beta),[1 1 nS]) ;
        % Combined spectral matrix
            Ft = sum(F.*Beta,3) ;
        % Joint transfer matrix
            [T,~] = eig(Ft) ;
        % Pole matrix
            PI = zeros([R nS]) ;
            for s = 1:nS
                PI(:,s) = diag(T\F(:,:,s)*T) ;
            end
        end
    % Signal Poles
        Z = reshape(PI,[R nS]).' ;
    % Wavevectors in the SHIFT basis
        shiftsCOS = any(this.isCOS.*shifts,2) ;
        K = zeros(size(Z)) ;
        K(~shiftsCOS,:) = log(Z(~shiftsCOS,:))/1i ; % FUNC = 'EXP' ;
        K(shiftsCOS,:) = acos(Z(shiftsCOS,:)) ; % FUNC = 'COS' ;
    % Wavevectors in the cartesian basis
        K = (shifts.*this.DECIM_K)\K ;
    end
end
    
    
%% AMPLITUDES/MODES ESTIMATION
methods
    function V = vandermondeMatrix(this,K,L) 
    % Build the vandermonde matrix of the signal
    % Here, decimation factors are no longer taken into account
        if nargin<2 ; K = this.extractPoles ; end
        if nargin<3 ; L = this.Lk ; end
        nD = numel(L) ;
    % Indices
        indV = zeros(prod(L),nD) ;
        for d = 1:nD
            indV(:,d) = repelem(repmat(0:L(d)-1,[1 prod(L(d+1:end))]),prod(L(1:d-1))) ;
        end
    % Wavevectors
        if any(this.isCOS) % Cosinuses searched, extended Vandermonde matrix
            K = kron( [1 (0.5-this.isCOS.')*2] ,K)  ; % [nD 2*nK]
        end
    % Vandermonde
        K = permute(K,[3 2 1]) ; % [1 nK nD]
        indV = permute(indV,[1 3 2]) ; % [prod(L) 1 nD]
        V = exp(1i*sum(indV.*K,3)) ; % [prod(L) nK]
    end
    
    function [A,v0] = extractAmplitudes(this,V)
    % Extract the exponential amplitudes from the signal
    % Here, decimation factors are no longer taken into account
        if nargin<2 ; V = this.vandermondeMatrix() ; end
    % Shift for the conditionning
        v0 = max(abs(V),[],1) ;
        V = V*diag(1./v0) ;
    % Amplitude estimation
        A = V\this.reshapedSignal ;
    % Vandermonde Shift compensation
        A = diag(1./v0)*A ;
    end
    
    function [U,PHI] = extractModes(this,A) 
    % Extract the modes from the signal
    % Here, decimation factors are no longer taken into account
        if nargin<2 ; A = this.extractAmplitudes() ; end
    % Reshaping
        lp = this.Lp ;
        nK = size(A,1) ;
        if ~any(this.isCOS) % No cosinuses, no phases to estimate
            U = reshape(A.',[lp nK]) ;
            if nargout>1 ; PHI = zeros(size(U)) ; end
        else
            Bplus = A(1:end/2,:) ;
            Bminus = A(end/2+1:end,:) ;
            U = reshape(2*sqrt(Bplus.*Bminus).',[lp nK]) ;
            PHI = reshape(1i/2*log(Bminus./Bplus).',[lp nK]) ;
        end
    end
end


%% SIGNAL MODEL
methods
    function SM = signalModel(this,V,A,reshp)
    % Compute the signal model from the extracted data
        if nargin<2 ; V = this.vandermondeMatrix ; end
        if nargin<3 ; A = this.extractAmplitudes(V) ; end 
        if nargin<4 ; reshp = true ; end
    % Model of this.reshapedSignal
        SM = V*A ;
    % Model of this.Signal
        if reshp
            SM = permute(reshape(SM,[this.Lk this.Lp]),this.permDims) ; 
        end
    end
end


%% UNCERTAINTY QUANTIFICATION
methods
    function opts = uncertaintyOptions(this)
    % Options for the estimation of uncertainties
    % HARD-CODED FOR NOW
        opts.hugeData = true ; % handles huge Data by dividing the data at points
        opts.estimate =   'std' ... % Standard Deviation
                        ... 'delta' ... % Sensibility or perturbation
                        ;
        opts.lin_method =     'conv' ... % Linearization method for the variance:
                            ... 'kron'...
                            ... 'none'...  
                            ;
        opts.covar_estim =   ... 'uniform' ... % Estimation of the perturbation covariance: 
                             'diagonal' ...
                            ;
        opts.formulation =   ... 'analytic' ... % with Vandermonde matrices etc. Do not work with cosinuses
                             'eigenspace' ... % from the eigenspace
                            ;
        opts.debug = this.DEBUG ;
    end
    
    function dS = modelError(this,V,A)
    % The signal error (or perturbation)
        if nargin<2 ; V = vandermondeMatrix(this) ; end
        if nargin<3 ; A = extractAmplitudes(this,V) ; end
        dS = this.reshapedSignal - this.signalModel(V,A,false) ;
    end
    
    function var_dS = errorVariance(this,dS,estim)
    % Estimate a covariance of the signal error (or perturbation)
        if nargin<2 ; dS = this.modelError() ; end
        if nargin<3 ; estim = this.uncertaintyOptions().covar_estim ; end
    % Standard deviation estimation: std
        switch estim
            case 'uniform'
                var_dS = var(dS(:),0); % scalar
            case 'diagonal'
                var_dS = abs(dS).^2 ; % matrix !
        end
    end
    
    function dK = poleUncertainty(this)
    % Uncertainty associated to the estimation of wavevectors
    end
    
    function dU = modeUncertainty(this)
    % Uncertainty associated to the estimation of amplitudes/modes
    end
    
    function [dK,dU] = computeUncertainties(this,K,V,A)
    % Return the uncertainties associated to each parameter
        if nargin<2 ; K = this.extractPoles ; end
        if nargin<3 ; V = this.vandermondeMatrix(this.K) ; end
        if nargin<4 ; A = this.extractAmplitudes(V) ; end
        dK = [] ; dU = [] ;
    % OPTIONS (hard-coded for now)
        opts = this.uncertaintyOptions() ;
    % COMMON DATA
        k = this.Kk ; m = this.Mk ; l = this.Lk ;
        decim = this.DECIM_K ;
        iHbH = this.blockHankelIndices() ;
        indP = this.modeIndices ;
        [nD,R] = size(K) ;
        if ismember(opts.formulation,{'eigenspace'}) 
            [W,lambda] = this.signalSubspace(R) ; 
            [~,T] = extractPoles(this,W,R) ;
        end
    % SIGNAL PERTURBATION AND COVARIANCE
        dS = this.modelError(V,A,false) ;
        if ~strcmp(opts.lin_method,'none') && strcmp(opts.estimate,'std')
            if(opts.debug) ; disp('        data covariance') ; end
            var_dS = errorVariance(this,dS,opts.covar_estim) ;
        end
    % WAVEVECTOR PERTURBATION
        if this.COMPUTE_dK
            if(opts.debug) ; disp('        wavevector uncertainties') ; end
            % Perturbed HbH matrix if needed
                if strcmp(opts.lin_method,'none') %|| strcmp(estimate,'delta')
                    if(opts.debug) ; disp('             dHss') ; end
                    dH = this.signalMatrix(dS) ;
                end
            % Selection matrices if needed
                if(opts.debug) ; disp('             selection matrices') ; end
                switch opts.lin_method
                    case 'kron' % M Matrix
                        Ip = speye(length(indP)) ;
                        I = speye(prod(l)) ;
                        M = sparse(prod(k)*prod(m)*length(indP),prod(l)*length(indP)) ;
                        for m = 1:prod(m)
                            Im = I(iHbH{1}(:,m),:) ;
                            for i = 2:numel(iHbH)
                                Im = Im + I(iHbH{i}(:,m),:) ;
                            end
                            Mm = kron(Ip,Im);
                            M((1:prod(k)*length(indP))+(m-1)*(prod(k)*length(indP)),:) = Mm ;
                        end
                    otherwise % Jdelta matrix (decimation selection matrix)
                        I = speye(prod(l),prod(l-m+1)) ;
                        Jdelta = I(iHbH{1}(:,1),:) ;
                        for i = 2:numel(iHbH)
                            Jdelta = Jdelta/2 + I(iHbH{i}(:,1),:)/2 ;
                        end
                end
            % INITIALIZATION
                dK = zeros(nD,R) ;
                shifts = eye(nD) ; % /!\ SHIFTS IS DEFAULT HERE !
            % HUGE DATA LOOP
                if opts.hugeData 
                    indP_bkp = indP(:)' ; % the for-loop will iterate
                else
                    indP_bkp = indP(:) ; % the for-loop will NOT iterate
                end
            % GOOOOOO !
                for indP = indP_bkp
                % Pre-Computations
                    % Right vector of the bilinear form
                        if(opts.debug) ; disp('             right eigenvectors of Css') ; end
                        switch opts.formulation
                            case 'analytic'
                            % Partial Vandermonde matrices
                                P = V(iHbH{1}(:,1),:) ;
                                Q = V(iHbH{1}(1,:),:) ;
                            % Complete right-Vandermonde Matrix
                                QA = zeros(prod(m)*length(indP),R) ;
                                for p = 1:length(indP) 
                                    QA((1:prod(m))+(p-1)*prod(m),:) = Q*diag(A(:,indP(p))) ;
                                end
                                x = (QA\eye(size(QA,1)))' ;
                            case 'eigenspace'
                                Hs = this.signalMatrix(this.signalModel(V,A,false),iHbH) ;
                                x = conj(((Hs'*W)./lambda(1:R).')*T) ;
                        end
                % Loop over the shifts
                    if(opts.debug) ; disp('             uncertainties') ; end
                    for s = 1:nD
                        [Jup,Jdwn] = this.selectionMatrices(shifts(s,:)) ;
                    % pre-compute
                        switch opts.formulation
                            case 'analytic'
                                vn = (Jup*P)\speye(size(Jdwn,1)) ;
                            case 'eigenspace'
                                vn = T\((Jup*W)\speye(size(Jdwn,1))) ;
                        end
                    % Loop over the R components
                        for r = 1:R
                            % The polar component
                                switch this.FUNC{s}
                                    case 'exp'
                                        PIrn = exp(1i*K(s,r)*decim(s)) ;
                                        Arn = exp(1i*K(s,r)*decim(s)) ;
                                    case 'cos'
                                        PIrn = cos(K(s,r)*decim(s)) ;
                                        Arn = abs(sin(K(s,r)*decim(s))) ;
                                end
                            % Linearization method
                                switch opts.lin_method
                                    case 'none' % NO LINEARIZATION, uncertainty only
                                        vrn = (vn(r,:)*(Jdwn-PIrn*Jup))' ;
                                        dK(s,r) = dK(s,r) + abs(vrn.'*dH*x(:,r)) ;
                                    case 'conv' % LINEAR / BY CONVOLUTION (USE OF THE HANKEL SHAPE OF Hss)
                                        vrn = (vn(r,:)*(Jdwn-PIrn*Jup)*Jdelta)' ;
                                        VRN = full(reshape(vrn,[1 l-m+1])) ;
                                        XR = reshape(x(:,r),[length(indP) m]) ;
                                        ZN = ifftn(bsxfun(@times,fftn(VRN,[1 l]),fftn(XR,[length(indP) l]))) ; % ND convolution
                                        zn = conj(ZN(:)) ;
                                    case 'kron' % LINEAR / BY VECTORIZATION (USES vec(A*X*B) = kron(B.',A)*vec(X) ) 
                                        vrn = (vn(r,:)*(Jdwn-PIrn*Jup))' ;
                                        zn = (kron((x(:,r)),vrn)'*M).' ;
                                end
                            % Perturbation estimate
                                if ~strcmp(opts.lin_method,'none')
                                    switch opts.estimate
                                        case 'std'
                                            switch opts.covar_estim
                                                case 'uniform'
                                                    dK(s,r) = dK(s,r) + var_dS*sum(abs(zn).^2) ;
                                                case 'diagonal'
                                                    dK(s,r) = dK(s,r) + sum(abs(zn).^2.*reshape(var_dS(indP,:),[],1)) ;
                                            end
                                        case 'delta'
                                            dK(s,r) = dK(s,r) + abs(zn)'*reshape(dS(indP,:),[],1) ;
                                    end
                                end
                        end % end of this signal order r
                    end % end of this shift s
                end % end of the points indP
            % Common terms
                if strcmp(opts.estimate,'std') ; dK = sqrt(dK) ; end
                dK = diag(decim)*abs(dK)/abs(Arn) ;
        end
    % AMPLITUDES PERTURBATION
        if this.COMPUTE_dU
            [u,phi] = this.extractModes(A) ;
            if nD>1 || any(this.isCOS(:)) ; dU = NaN*ones(size(u)) ; return ; end
            if(opts.debug) ; disp('        amplitudes uncertainty') ; end
            % Pre-compute the inverse of V
                if(opts.debug) ; disp('             inverse of V') ; end
                invV = V\eye(prod(l)) ;
            % Uncertainty
                switch opts.estimate
                    case 'std'
                        dU = zeros(R,prod(this.Lp));
                        for r = 1:R
                            zrn = invV(r,:)' ;
                                switch opts.covar_estim
                                    case 'uniform'
                                        dU(r,:) = sqrt(var_dS*sum(abs(zrn).^2)) ;
                                    case 'diagonal'
                                        dU(r,:) = sqrt(var_dS*abs(zrn(:)).^2).' ;
                                end
                            dU = dU*2 ; % Yes, I don't know why this factor 2...
                        end
                    case 'delta'
                        dV = ((0:l-1)'.*V)*diag(dK(1,:)) ;
                        dU = abs(abs(invV)*(dS.' + abs(dV)*u.'));
                end
            % Final Processing
                dU = reshape(dU.',[this.Lp R]) ;
        end
    end
end
    
    
end





