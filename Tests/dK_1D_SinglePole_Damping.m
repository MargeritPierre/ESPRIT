%% SINGLE POLE CASE , 1D, INFLUENCE OF THE DAMPING

clc
clear all
clf 

L = 100 ;
K = 1:L-1 ;
M = L-K+1 ;
 
omega = .02 ;
alpha = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.14 0.16] ; linspace(0,.2,11) ;logspace(-3,-1,5) ;
txt = text() ; txt(1) = [] ;
pl = plot(0,0) ; pl(1) = [] ;
minPts = plot(0,0) ; minPts(1) = [] ;
labelPos = round(L/4) ;
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
        y = v(1:k-1)\(Jdown-z*Jup) ; 
        e = conv(x,y) ;
        Err(i) = norm(e)^2 ;
    end
    data = 1./(Err/CRB/2/abs(z)^2) ; data(isinf(data))=NaN ;
    pl(a) = plot(M/L,data) ;
    txt(a) = text(pl(a).XData(labelPos),pl(a).YData(labelPos),[num2str(alpha(a)),'']);
    [~,i] = max(data(1:floor(end/2))) ;
    minPts(a) = plot(M(i)/L,data(i),'.k') ;
    set(gca,'yscale','log')
    grid on
end
set(pl,'color','k','linewidth',1.0)
set(minPts,'markersize',15)
set(txt,'backgroundcolor','w','fontsize',18)
set(gca,'ylim',sort(1./[1 100]))

xlabel '\boldmath$r$'
%ylabel '\boldmath$\eta$'
title ' \boldmath $\mathrm{CRB}_\omega/\mathrm{var}(\omega)$'




%% ANALYTICAL VALUES IN THE CONSERVATIVE CASE

clc
clear all
clf 

L = 100 ;
K = 2:L-1 ;
M = L-K+1 ;

z = exp(0.05i) ;

kk = linspace(min(K),max(K),1000) ;
mm = L-kk+1 ;
E = 2.*min(kk-1,mm)./(kk-1).^2./mm.^2 ;

CRB = 12/L^3 ;

plot(kk/L,E/CRB,'k','linewidth',1) ;
set(gca,'yscale','log')
grid on

Err = zeros(length(K),1) ;
v = z.^(0:L-1).' ;
for i = 1:length(K)
    k = K(i) ; m = M(i) ;
    x = conj(v(1:m))/norm(v(1:m))^2 ;
    Jup = [eye(k-1) zeros(k-1,1)] ;
    Jdown = [zeros(k-1,1) eye(k-1)] ;
    y = zeros(1,k) ;
    y(1) = -z ;
    y(end) = conj(z)^(k-2) ;
    y(2:end-1) = conj(z).^(0:k-3).'*(1-abs(z)^2) ;
    y = y/norm(v(1:k-1))^2 ;
    %y = v(1:k-1)'*(Jdown-z*Jup)/norm(v(1:k-1))^2 ; %(v(2:k+1)'-z*v(1:k)').'/norm(v(1:k-1))^2 ; %(Jup*v(1:k))\(Jdown-z*Jup) ; %
    e = conv(x,y) ;
    Err(i) = norm(e)^2 ;
end

plot(K/L,Err/CRB,'or','markersize',3,'linewidth',1.0)
