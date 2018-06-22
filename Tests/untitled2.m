

clc
clear all
%clf 

L = 100 ;

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 1 ; ...
           'R0' , 2 ; ...
           'DECIM' , [1 1] ; ...
           'FUNC' , 'exp' ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'K0' , [] ; ...
          }' ;
f = 0.05 ;
df = linspace(-1,1,100)*.1 ;
z = exp(1i*f) ;
a = [1 1] ;
SNR = 1e3 ;
nMCMC = 10 ;
profiler = false ;

K = zeros(length(df),2) ;
dK = K ;
wtbr = waitbar(0) ;
ti = tic ;
if profiler ; profile on ; end
for s = 1:length(df)
    for m = 1:nMCMC
        dz = exp(1i*df(s)) ;
        V = [(z*dz).^(0:L-1).' (z*conj(dz)).^(0:L-1).'] ;
        Amp = ((rand(1,2)*2-1)+1i*(rand(1,2)*2-1)) ;
        Amp = Amp./abs(Amp).*a ;
        Signal = V*Amp(:) ;
        Rn = rand(size(Signal))*2-1 ;
        In = rand(size(Signal))*2-1 ;
        noise = (Rn+1i*In) ;
        noise = noise*norm(Signal)/norm(noise)/SNR ;
        out = ESPRIT(Signal+noise,EspArgs{:}) ;
        [~,ind] = sort(real(out.K)) ;
        K(s,:,m) = out.K(ind) ;
        dK(s,:,m) = out.dK(ind) ;
        if toc(ti)>.2
            ti = tic ;
            wtbr = waitbar((m+(s-1)*nMCMC)/nMCMC/length(df),wtbr) ;
        end
    end
end
delete(wtbr) ;
drawnow ;



varK = var(K,1,3) ;

plot(df,mean(dK,3),'k') ;
plot(df,varK,'.','markersize',15) ;

set(gca,'yscale','log')
grid on



%%
plot(df,Err,'k')
set(gca,'yscale','log')
grid on


%%

clc
clear all
clf 

L = 100 ;
M = round(2*L/3) ;
K = L-M+1 ;

f = 0.05 ;
df = linspace(-1,1,1000)*.1 ;
z = exp(1i*f) ;

Err = zeros(length(df),2) ;
for i = 1:length(df)
    dz = exp(1i*df(i)) ;
    v = [(z*dz).^(0:L-1).' (z*conj(dz)).^(0:L-1).'] ;
    x = conj(v(1:M,:))/norm(v(1:M,:))^2 ;
    Jup = [eye(K-1) zeros(K-1,1)] ;
    Jdown = [zeros(K-1,1) eye(K-1)] ; 
    for r = 1:2
        y = (v(1:K-1,:)'*(Jdown-z*exp(1i*df(i)*(1.5-r)*2)*Jup)/norm(v(1:K-1))^2).' ;
        e = conv(x(:,r),y(:,r)) ;
        Err(i,r) = norm(e)^2 ;
    end
end


plot(df,Err,'k')
set(gca,'yscale','log')
grid on


