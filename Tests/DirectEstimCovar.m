clc
clf ;

U = Signal(:,indT).' ;
Cuu = U'*U ;

[Wr,sig] = eig(Cuu,'vector') ;
[sig,ind] = sort(sig,'descend') ;
Wr = Wr(:,ind) ;

MDL = (sig(1:end-1)-sig(2:end))./sig(1:end-1) ;
plot(MDL)
%set(gca,'yscale','log') ;
axis tight

W = bsxfun(@(w,s)w./s, U*Wr, sqrt(sig(:)).') ;


%% ESPRIT-LIKE

R0 = 1:200 ;
clf
fftU = log10(mean(abs(fft(U,[],1)),2)) ;
fftU = (fftU-min(fftU))/(max(fftU)-min(fftU))*(max(R0)-min(R0))+min(R0);
fftFreq = (0:length(indT)-1)/length(indT)*EXP.Fe ;
plot(fftFreq,fftU) ;

set(gca,'xscale','log') ;
%set(gca,'yscale','log') ;
set(gca,'xlim',[fftFreq(2) EXP.Fe/2])

pts = plot(NaN*zeros(length(R0).',max(R0)),NaN*zeros(length(R0),max(R0)).','.k') ;

for rr = 1:length(R0) ;
    Wup = W(1:end-1,1:R0(rr)) ;
    Wdwn = W(2:end,1:R0(rr)) ;
    F = Wup\Wdwn ;
    f = -1i*log(eig(F))/2/pi/EXP.dt ;
    pts(rr).XData(1:R0(rr)) = abs(real(f)) ;
    pts(rr).YData(1:R0(rr)) = R0(rr) ;
    drawnow ; 
end


%% MODES
 clf ; 
 m = 1 ;
 trisurf(EXP.Tri,EXP.X0,EXP.Y0,real(Wr(:,1)))
 axis equal
 axis tight







