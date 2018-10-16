clc
clear all

L = 15 ;
d = 2 ;

K = 5 ;
M = L-(K-1)*d ;

s = 1:L ;

H = hankel(0:d*(K-1),d*(K-1):L-1)+1 ;
H = H(1:d:end,:)

x = H' ;
v = H*x

v2 = convn(s(:),flip(x),'full') ;
v2 = v2(M+(0:d:(K-1)*d),:) 
norm(v-v2)


