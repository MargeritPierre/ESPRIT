L = unique(round(logspace(1,4,50))) ;
T = .1 ;

tFFT = [] ;
tCONV = [] ;

wtbr = waitbar(0) ;
for l=1:length(L)
    % Generate Signals
        s1 = (rand(L(l),1)*2-1) + 1i*(rand(L(l),1)*2-1) ;
        s2 = (rand(L(l),1)*2-1) + 1i*(rand(L(l),1)*2-1) ;
    % CONV Version
        t = tic ;
        it = 0 ;
        while toc(t)<T
            Cconv = conv(s1,s2,'full') ;
            it=it+1 ;
        end
        tCONV(l) = toc(t)/it ;
    % FFT Version
        t = tic ;
        it = 0 ;
        while toc(t)<T
            Cfft = ifft(fft(s1,2*L(l)-1).*fft(s2,2*L(l)-1)) ;
            it=it+1 ;
        end
        tFFT(l) = toc(t)/it ;
    % Waitbar
        wtbr = waitbar(l/length(L),wtbr) ;
end
delete(wtbr) ;
%%
clf ; 
plot(L,tCONV,'r') ;
plot(L,tFFT,'b') ;
set(gca,'xscale','log','yscale','log')
grid on

% Linear regressions
fCONV = L.^2 ;
fFFT = L.*log(L) ;
aCONV = fCONV(:)\tCONV(:) 
aFFT = fFFT(:)\tFFT(:) 
plot(L,aCONV*fCONV,':r') ;
plot(L,aFFT*fFFT,':b') ;




