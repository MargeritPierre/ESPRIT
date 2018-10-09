%% ANALYTICAL RESULTS TEST, AS A FUNCTION OF THE SPATIAL SMOOTHING

clc
%clf
clear all

L = 1000 ;
d = 1  ;

z = exp(1i*(0.02+0.05i*100/L)) ;

M = ceil(L/2):L ;
K = (L-M)./d + 1 ;
K = unique(round(K)) ;
M = L - (K-1).*d ;

K(M==0)=[] ;
M(M==0)=[] ;

g = zeros(1,length(K)) ;
gf = zeros(1,length(K)) ;
for k = 1:length(K) ; length(K) ; 
    v1M = z.^(0:M(k)-1).' ;
    vdK = z.^(0:d:(K(k)-1)*d).' ;
    J = eye(K(k)*d) ;
    J = J(1:d:end,:) ;
    NORM = norm(v1M)^4*norm(vdK)^4 ;
    % Raw computation
        g(k) = norm(conv(v1M,J.'*vdK))^2/NORM ;
    % Step 1
        vdKd = zeros(K(k)*d,1) ;
        vdKd(1:d:end) = z.^(0:d:K(k)*d-1).' ;
        %gf(k) = norm(conv(v1M,vdKd))^2/NORM ;
    % Step 2
%         vv = zeros(M(k)+K(k)*d-1,1) ;
%         for m = 0:K(k)-1
%             v0 = [zeros(m*d,1);z^(m*d);zeros((K(k)-m)*d-1,1)] ;
%             vv = vv + conv(v1M,v0) ;
%         end
        %gf(k) = norm(vv)^2/NORM ;
    % Step 3
%         vv = zeros(M(k)+K(k)*d-1,1) ;
%         for m = 0:K(k)-1
%             v0 = [zeros(m*d,1);v1M;zeros((K(k)-m)*d-1,1)] ;% 
%             vv = vv + z^(m*d)*v0 ;
%         end
%         gf(k) = norm(vv)^2/NORM ;
    % Step 4
%         zd = z^d ;
%         Np = max(K(k),(M(k)/d)) ;
%         Nm = min(K(k),(M(k)/d)) ; 
%         vt = [(1:Nm).*zd.^(0:Nm-1),Nm*zd.^(Nm:Np-1),(Np+Nm-1-(Np:Np+Nm-2)).*zd.^(Np:Np+Nm-2)].' ;
%         vd = z.^(0:d-1).' ;
%         gf(k) = norm(vd)^2*norm(vt)^2/NORM ;
    % Final expression
        %gf(k) = d/Np*(1-(Nm^2-1)/(3*Np*Nm)) ;
        %gf(k) = d*Nm^2*(Np-1/3/Nm*(Nm^2-1))/NORM ;
    % with spatial smoothing factor
        s = M(k)/L ;
        %gf(k) = 1/s/L*(1 - (1-s)/3/s*((1-s)+2*d/L)/((1-s)+d/L)) ;
end

plot(M/L,g*L)
ax = gca ; ax.ColorOrderIndex=ax.ColorOrderIndex-1 ;
plot(M/L,gf*L,'o')
set(gca,'yscale','log')
axis tight
grid on


%% FIXED SPATIAL SMOOTHING, AS A FUNTION OF THE DAMPING

clc
clf
clear all

L = 100 ;
D = unique(ceil([1 5 10 20]/100*L)) 
s = 2/3 ;
alpha = -linspace(0,10,50) ;



M = s*L ;
K = (L-M)./D + 1 ;
K = round(K)
M = L - (K-1).*D 

g = zeros(length(alpha),length(D)) ;
for a = 1:length(alpha)
    for k = 1:length(K) ;
        z = exp(alpha(a)/L) ;
        d = D(k) ;
        v1M = z.^(0:M(k)-1).' ;
        vdK = z.^(0:d:(K(k)-1)*d).' ;
        J = eye(K(k)*d) ;
        J = J(1:d:end,:) ;
        NORM = norm(v1M)^4*norm(vdK)^4 ;
        g(a,k) = norm(conv(v1M,J.'*vdK))^2/NORM ;
    end
end

plot(alpha(:),g*L)
set(gca,'yscale','log')
%axis tight
grid on




