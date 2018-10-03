%% TEST THE PERFORMANCE OF THE ESTIMATION OF ONE EXPONENTIAL
% ANALYTICAL RESULTS

clc
clear all
%close all

L = 1000 ;
delta = 2 ;
k = 0.01+0.002i ;
nu = exp(1i*k) ;

K = floor((0:L)/delta) ; K = K(K>1) ;
M = L-(K-1).*delta ;
r = M./L ;

CRB = 6./L^3 ;

vand = @(delta,k)(nu^delta).^(0:k-1).' ;

f = zeros(1,length(r))*NaN ;
for rr = 1:length(r)
    if K(rr)==1 ; continue ; end
    v_1_M = vand(1,M(rr)) ;
    v_delta_Kminus1 = vand(delta,K(rr)-1) ;
    Jdelta = speye(delta*(K(rr)-1)+1) ; Jdelta = Jdelta(1:delta:end,:) ;
    Jup = speye(K(rr)) ; Jup = Jup(1:end-1,:) ;
    Jdwn = speye(K(rr)) ; Jdwn = Jdwn(2:end,:) ;
    y = (Jdelta'*(Jdwn'-conj(nu^delta)*Jup')*v_delta_Kminus1) ;
    zi = conv(v_1_M,full(y)) ;
    f(rr) = norm(zi)^2/norm(v_1_M)^4/norm(v_delta_Kminus1)^4 ;
end



%fig = clf ;
plot(r,f./CRB/2/delta^2,'k')
set(gca','yscale','log')
set(gca,'ylim',[1 100])
grid on



