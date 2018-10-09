%% SINGLE POLE CASE , 1D, INFLUENCE OF THE DAMPING

clc
clear all
%clf 

L = 200 ;
d = 20 ; 
linestyle = ':'

M = linspace(2,L,100) ;
K = (L-M)/d+1 ;
K = unique(round(K)) ;
M = L-(K-1)*d ;
K(M==0)=[] ;
M(M==0)=[] ;

 
omega = .03 ;
alpha = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1]*100/L ; linspace(0,.2,11) ;logspace(-3,-1,5) ;
txt = text() ; txt(1) = [] ;
pl = plot(0,0) ; pl(1) = [] ;
minPts = plot(0,0) ; minPts(1) = [] ;
labelPos = round(length(K)/5) ;
for a = 1:length(alpha) ;
    z = 1*exp(1i*omega+alpha(a)) ;
    Err = zeros(length(K),1) ;
    v = z.^(0:L-1).' ;
    Pwr = norm(v)^2/L ;
    CRB = 6/L^3/Pwr ;
    for i = 1:length(K)
        k = K(i) ; m = M(i) ;
        x = conj(v(1:m))/norm(v(1:m))^2 ;
        Jup = [eye(k-1) zeros(k-1,1)] ;
        Jdown = [zeros(k-1,1) eye(k-1)] ;
        Jdelta = eye(k*d) ; Jdelta = Jdelta(1:d:end,:) ;
        y = v(1:d:(k-1)*d)\(Jdown-z^d*Jup)*Jdelta ; 
        e = conv(x,y) ;
        Err(i) = norm(e)^2 ;
    end
    data = 1/d^2.*Err/CRB/2/abs(z^d)^2 ; data(isinf(data))=NaN ;
    pl(a) = plot(M/L,data,linestyle) ;
    txt(a) = text(pl(a).XData(labelPos),pl(a).YData(labelPos),[num2str(alpha(a)*L),'']);
    [~,i] = max(data(1:floor(end/2))) ;
    %minPts(a) = plot(M(i)/L,data(i),'.k') ;
    set(gca,'yscale','log')
    grid on
end
set(pl,'color','k','linewidth',1.0)
set(minPts,'markersize',15)
set(txt,'backgroundcolor','w','fontsize',18)
set(gca,'ylim',sort([1 100]))

xlabel '\boldmath$r$'
%ylabel '\boldmath$\eta$'
title ' \boldmath $\mathrm{CRB}_\omega/\mathrm{var}(\omega)$'

delete(minPts)


%% SINGLE POLE CASE , 1D, MULTIPLE DECIMATION,MULTIPLE DAMPING, 
%fixed spatial smooting factor s

clc
clear all
clf 

L = 1000 ;
D = max(round([1 2 5 10 15 50 100 200]/1000*L),1) ;
alpha = linspace(0,10,150)/L ;
s = 2/3 ;
linestyle = '-' ;

M = round(s*L) ;
K = (L-M)./D+1 ;
K = floor(K) 
M = L-(K-1).*D 

 
omega = .03 ;
pl = plot(0,0) ; pl(1) = [] ;
DATA = zeros(length(alpha),length(D)) ;
for a = 1:length(alpha) ;
    z = 1*exp(1i*omega+alpha(a)) ;
    Err = zeros(length(K),1) ;
    v = z.^(0:L-1).' ;
    Pwr = norm(v)^2/L ;
    CRB = 6/L^3/Pwr ;
    for i = 1:length(D)
        k = K(i) ; m = M(i) ; d = D(i) ;
        x = conj(v(1:m))/norm(v(1:m))^2 ;
        Jup = [eye(k-1) zeros(k-1,1)] ;
        Jdown = [zeros(k-1,1) eye(k-1)] ;
        Jdelta = eye(k*d) ; Jdelta = Jdelta(1:d:end,:) ;
        y = v(1:d:(k-1)*d)\(Jdown-z^d*Jup)*Jdelta ; 
        e = conv(x,y) ;
        Err(i) = norm(e)^2 ;
        DATA(a,i) = 1/d^2.*Err(i)/CRB/2/abs(z^d)^2 ;
    end
end
pl = plot(alpha*L,DATA,linestyle) ;
%set(pl,'color','k','linewidth',1.0)
%set(gca,'ylim',sort([1 100]))
set(gca,'yscale','log')
grid on

%xlabel '\boldmath$r$'
%ylabel '\boldmath$\eta$'
%title ' \boldmath $\mathrm{CRB}_\omega/\mathrm{var}(\omega)$'







