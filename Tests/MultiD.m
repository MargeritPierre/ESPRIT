clc
clear all
close all

R = 5 ; 
D = 3 ;
Dc = [] ; [1 2 3] ;
De = setdiff(1:D,Dc) ;

L = repmat(15,[1 D]) ;

K = (rand(R,D)*2-1)*pi/2*.99+1i*(rand(R,D)*2-1)*pi/2*.1 ;
K(:,Dc) = acos(cos(K(:,Dc))) ;
U = (rand(R,1)*2-1) ;

FUNC = repmat('exp',[D 1]) ;
FUNC(Dc,:) = repmat('cos',[length(Dc) 1]) ;

X = [] ;
for d = 1:D
    xd = reshape(0:L(d)-1,[ones(1,d-1) L(d) ones(1,D-d) 1]) ;
    X = cat(D+1,X,repmat(xd,[L(1:d-1) 1 L(d+1:end)])) ;
end

XX = repmat(X,[ones(1,D+1) R]) ;
KK = repmat(reshape(K.',[ones(1,D) D R]),[L 1]) ;

isEXP = zeros(1,D) ; isEXP(De) = 1 ;
isCOS = zeros(1,D) ; isCOS(Dc) = 1 ;
isEXP = repmat(reshape(isEXP,[ones(1,D) D 1]),[L 1 R]) ;
isCOS = repmat(reshape(isCOS,[ones(1,D) D 1]),[L 1 R]) ;

EXP = exp(sum(1i.*XX.*KK.*isEXP,D+1)) ;
COS = cos(sum(XX.*KK.*isCOS,D+1)) ;
UU = repmat(reshape(U,[ones(1,D) 1 R]),[L 1]) ;

SIGNAL = reshape(sum(UU.*EXP.*COS,D+2),[L 1]) ;

EspArgs = {... Arguments for ESPRIT
           'DIMS_K' , 1:D ; ...
           'R0' , R ; ...
           'DECIM' , [ones(1,D) 1] ; ...
           'FUNC' , FUNC ; ...
           'FIT' , 'LS' ; ...
           'SOLVER' , 'eig' ; ...
           'DEBUG' , false ; ...
           'M/L' , [1/3] ; ...
          }' ;
t = tic ;      
out = ESPRIT(SIGNAL,EspArgs{:}) ;
toc(t) ;
Kesp = out.K.' ;

[~,indESP] = sort(real(Kesp(:,1))) ;
[~,ind] = sort(real(K(:,1))) ;

Kesp(indESP,:),K(ind,:)
mean(mean(abs(Kesp(indESP,:)-K(ind,:))))
out




      
  