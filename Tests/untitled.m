

clc
clear all
clf 

L = 100 ;
K = 1:L-1 ;

M = L-K+1 ;

z = exp(0.05i) ;

E = 2.*min(K-1,M)./(K-1).^2./M.^2 ;

plot(K/L,E,'k') ;
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

plot(K/L,Err,'or')


