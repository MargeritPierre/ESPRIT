clc
clear all
%close all


L = [1 1]*20 ;
fun = {'cos' 'exp'} ;
K = [1*(1+1i) 5].' ;
A = [1] ;

D = length(L) ;
de = ismember(fun,{'exp'}) ;
dc = ismember(fun,{'cos'}) ;
De = sum(de) ;
Dc = sum(dc) ;
R = size(K,2) ;

K = K./repmat(L',[1 R])*pi ;

X = [] ;
for d = 1:D
    X = cat(D+1,X,repmat(reshape(0:L(d)-1,[ones(1,d-1) L(d) ones(1,D-d)]),[L(1:d-1) 1 L(d+1:D)])) ;
end
X = permute(X,[2 1 (3:D+1)]) ;
Ke = reshape(diag(de)*K,[ones(1,D) D R]) ;
Kc = reshape(diag(dc)*K,[ones(1,D) D R]) ;
Ue = exp(1i*sum(bsxfun(@times,X,Ke),D+1)) ;
Uc = cos(sum(bsxfun(@times,X,Kc),D+1)) ;
U = sum(bsxfun(@times,reshape(A,[ones(1,D+1) R]),Ue.*Uc),D+2) ; 

fig = clf ;
fig.WindowStyle = 'docked' ;
colormap(linspecer(1000)*.9+.1);

srf = surf(X(:,:,2),X(:,:,1),real(U)...
            ,'facecolor','interp'...
            ,'linestyle','-.'...
            ) ;
    %shading interp ;
    axis tight ;
    myaxisequal('xy') ;
    set(gca,'view',[36 45])
    %axis off
    
    
    
    
    
    
%%    
    matlab2tikz('essai.tikz','checkForUpdates',false)
    
    
    
    
    