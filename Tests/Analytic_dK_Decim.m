clc
clear all
clf ;

L = 100 ;
decim = 1:10 ;
F = [0.01-0.0i] ;  [-0.05 0.04] ;

U = ones(size(F)) ;
t = 0:L-1 ;
Signal = U(:).'*exp(1i*F(:)*t(:)') ;

CRB = 12/L^3/(norm(Signal)^2/L) ;

lines = plot(0,0) ; lines(1) = [] ;
pts = plot(0,0) ; pts(1) = [] ;
for d = 1:length(decim)
    Q = decim(d) ;
    K = 2+length(F):floor((L-length(F))/Q)+1 ;
    M = zeros(length(K),1) ;
    varK = zeros(length(K),length(F)) ;
    for i = 1:length(K)
        EspArgs = {... Arguments for ESPRIT
                   'DIMS_K' , 2 ; ...
                   'R0' , length(F) ; ...
                   'DECIM' , [1 Q] ; ...
                   'FUNC' , 'exp' ; ...
                   'FIT' , 'LS' ; ...
                   'SOLVER' , 'eig' ; ...
                   'DEBUG' , false ; ...
                   'K0' , [K(i)] ; ...
                  }' ;
        out = ESPRIT(Signal,EspArgs{:}) ;
        [~,ind] = sort(real(out.K)) ;
        varK(i,:) = out.dK(ind) ;
        M(i) = out.Mk ;
    end
    pts = [pts plot(M,varK/CRB)] ;
    if 0
        kk = linspace(K(1),K(end),1000) ;
        mm = L-(kk-1)*Q ;
        Err = 2*min(mm,(kk-1)*Q)./(((kk-1)*Q).^2.*mm.^2) ;
        lines(d) = plot(kk/L,Err/CRB) ;
    end
end
set(lines,'linewidth',1) ;
%set(pts,'linestyle',':','marker','.','markersize',15,'linewidth',1) ;
set(pts,'linewidth',1) ;

%set(gca,'xscale','log')
set(gca,'yscale','log')
grid on
