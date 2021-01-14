%% DEMONSTRATION SCRIPTS FOR THE USE OF THE ESPRIT FUNCTION

%% 1D SERIES OF EXPONENTIALS: S(t) = sum_r { A_r * exp(1i*w_r*t) }
    clearvars
    % Parameters
        nT = 100 ; % signal length
        dt = 0.1 ; % time step
        R = 5 ; % exact signal order < nT/2
        maxFreq = 0.25*pi/dt ; % maximum absolute frequency < pi/dt
        maxDamp = 0.05 ; % maximum damping << 1
        R0 = ... R ... known signal order
             1:min(2*R,floor(nT/2)) ... signal order candidates
             ;
        snr = 1e1 ; % signal to noise ratio
    % Signal construction
        t = (0:nT-1)*dt ; % time
        w = 2*(rand(1,R)-0.5)*maxFreq ; % random pulsations
        w = w.*(1+2i*(rand(1,R)-0.5)*maxDamp) ; % random damping
        A = rand(1,R).*exp(2i*pi*rand(1,R)) ;  % random complex amplitudes
        S = sum(A.*exp(1i*w.*t(:)),2) ; % signal [nT 1]
    % Additive noise
        b = rand(nT,1) ;
        S = S + b*(norm(S)/norm(b)/snr) ;
    % Apply ESPRIT
        esp = ESPRIT(S,'R0',R0,'COMPUTE_U',true) ; % amplitudes are not extracted by default
    % Signal model
        tt = linspace(t(1),t(end),10000) ;
        Se = sum(esp.U(:).'.*exp(1i*esp.K(:).'.*tt(:)/dt),2) ;
    % Compare results
        if numel(w)==numel(esp.K) % if signal orders coincide, compare
            [~,iw] = sort(real(w)) ;
            [~,iwe] = sort(real(esp.K)) ;
            errW = abs(esp.K(iwe)/dt - w(iw)) 
            errA = abs(esp.U(iwe) - A(iw))
        end
    % Compare signal and signal model
        clf ; hold on ; 
        plot3(real(S),imag(S),t,'.') ; % samples
        plot3(real(Se),imag(Se),tt,'-','linewidth',1) ; % signal model
        xlabel 'real(S)' ; ylabel 'imag(S)' ; zlabel 'time'
        view(45,15) ;
        
        
%% 2D SERIES OF EXPONENTIALS: S(x) = sum_r { A_r * exp(1i*k_r.x) }
    % Parameters
        L = 20*[2 1.5] ; % signal length
        dx = [0.1 0.2] ; % time step
        R = 10 ; % exact signal order < prod(L/2)
        R0 =  R ... known signal order
             ... 1:min(2*R,floor(prod(L/2))) ... signal order candidates
             ;
        snr = 5e2 ; % signal to noise ratio
        maxKi = 0.25*pi./dx ; % maximum absolute wavevector < pi./dx
        maxDamp = 0.05 ; % maximum damping << 1
    % Uniform grid of "measurement" points
        X = arrayfun(@colon,L*0,L-1,'UniformOutput',false) ;
        [X{:}] = ndgrid(X{:}) ;
        X = cat(numel(L)+1,X{:}) ;
        X = reshape(X,[],numel(L)).*dx ;
    % Signal construction
        k = 2*(rand(2,R)-0.5).*maxKi(:) ;
        k = k.*(1 + 2i*(rand(2,R)-0.5)*maxDamp(:)) ;
        A = rand(1,R).*exp(2i*pi*rand(1,R)) ;  % random complex amplitudes
        S = sum(A.*exp(1i*X*k),2) ; % signal [prod(L) 1]
        S = reshape(S,L) ;
    % Additive noise
        b = rand(L) ;
        S = S + b*(norm(S(:))/norm(b(:))/snr) ;
    % Apply ESPRIT
        esp = ESPRIT(S,[1 2],R0,'COMPUTE_U',true) ; % amplitudes are not extracted by default
    % Signal model
        Se = sum(esp.U.*exp(1i*X*(esp.K./dx(:))),2) ; % signal [prod(L) 1]
    % Compare results
        if numel(k)==numel(esp.K) % if signal orders coincide, compare
            [~,ik] = sort(real(k(1,:))) ;
            [~,ike] = sort(real(esp.K(1,:))) ;
            errK = abs(esp.K(:,ike).'./dx - k(:,ik).') 
            errA = abs(esp.U(ike) - A(ik))
        end
    % Compare signal and signal model
        clf ; hold on ; 
        plot3(X(:,1),X(:,2),real(S(:)),'.k') ; % samples
        surf(reshape(X(:,1),L),reshape(X(:,2),L),reshape(real(Se),L),'facecolor','interp','edgecolor',[1 1 1]*0.5) ; % signal model
        xlabel 'x' ; ylabel 'y' ; zlabel 'S'
        view(45,60) ;







    