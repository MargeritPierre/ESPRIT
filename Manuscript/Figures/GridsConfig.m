%% SUB GRIDS SAMPLES FOR THE 2D CASE WITH EXPONENTIAL ONLY

clc
clear all
close all ;

L = [8 6] ;

Q = [1 0 ; 0 1 ; 3 3 ; -2 1].' ;

ii = repmat((0:L(1)-1)',[1 L(2)]) ;
jj = repmat((0:L(2)-1),[L(1) 1]) ;

colors = linspecer(size(Q,2))*.9 ;
markersize = 15/1.6 ;
linewidth = 1.5 ;
spacing = .08 ;
l = norm(L) ;

fig = clf ;
fig.Position(3:4) = [1029 198] ;

nQ = size(Q,2) ;
for q = 1:nQ
    Dwn = (ii-Q(1,q)<L(1)) & (ii-Q(1,q)>-1) & (jj-Q(2,q)<L(2)) & (jj-Q(2,q)>-1) ;
    Up = (ii+Q(1,q)<L(1)) & (ii+Q(1,q)>-1) & (jj+Q(2,q)<L(2)) & (jj+Q(2,q)>-1) ;
    ax = mysubplot(1,nQ,q) ;
        plot(ii(:),jj(:),'.','markersize',5,'linewidth',1,'color',[1 1 1]*0.)
        plot(ii(Dwn),jj(Dwn),'o','color',colors(q,:),'markersize',markersize,'linewidth',linewidth) ;
        plot(ii(Up),jj(Up),'.','color',colors(q,:),'markersize',markersize*1.2,'linewidth',linewidth) ;
        axis equal
        axis tight
        ax.XLim = [0 L(1)-1] + spacing*[-1 1]*l ;
        ax.YLim = [0 L(2)-1] + 1*spacing*[-1 1]*l ;
        ax.XTick = [] ; ax.YTick = [] ;
        axis off
        xlabel(['(',char(96+q),')', ' $\vect q = [',num2str(Q(1,q)),'\;\;',num2str(Q(2,q)),']$' ],'visible','on') ;
        %ax.OuterPosition = [(q-1)*1/nQ 0 1/nQ 1] ;
end

% SAVE THE TIKZ FIGURE 
    cleanfigure ;
    figPos = get(gcf,'position') ;
    normHeight = figPos(4)/figPos(3) ;
    matlab2tikz('SubGridsEXP.tikz'...
                ,'width','\figwidth'...
                ,'height',[num2str(normHeight,4),'\figwidth']...
                ...,'noSize',true...
                ...,'parseStringsAsMath',true...
                ) ;
            
            
%% SUB GRIDS SAMPLE FOR THE 2D CASE WITH EXP AND COS

clc
clear all
close all ;

L = [8 6] ;

Q = [1 0 ; 0 1 ; -2 3 ; 1 1 ; 1 2 ; -1 2].' ;
fun = {['exp';'cos'],['exp';'cos'],['exp';'exp'],['cos';'cos'],['cos';'cos'],['cos';'cos']} ;


D = length(L) ;
nQ = size(Q,2) ;

ii = [] ;
for d = 1:D
    ii = cat(D+1,ii,repmat(reshape(0:L(d)-1,[ones(1,d-1) L(d) ones(1,D-d)]),[L(1:d-1) 1 L(d+1:end)])) ;
end

colors = linspecer(size(Q,2))*.9 ;
marker = 'o' ;
markersize = 7 ; 14 ;
linewidth = 1.5 ;
spacing = .08 ;
l = norm(L) ;

for q = 1:nQ
fig = figure ;
fig.Position(3:4) = [1081 188] ;

ind = permute(ii,[D+1 1:D]) ;
Up = true(1,prod(L)) ;
Dwn1 = true(1,prod(L)) ;
Dwn2 = true(1,prod(L)) ;
for d = 1:D
    switch fun{q}(d,:)
        case 'exp'
            DwnD = (ind(d,:)-Q(d,q)<L(d)) & (ind(d,:)-Q(d,q)>-1) ;
            Dwn1 = Dwn1 & DwnD ;
            Dwn2 = Dwn2 & DwnD ;
            Up = Up & (ind(d,:)+Q(d,q)<L(d)) & (ind(d,:)+Q(d,q)>-1) ;
        case 'cos'
            Dwn1 = Dwn1 & (ind(d,:)-2*Q(d,q)<L(d)) & (ind(d,:)-2*Q(d,q)>-1) ;
            Dwn2 = Dwn2 & (ind(d,:)+2*Q(d,q)<L(d)) & (ind(d,:)+2*Q(d,q)>-1) ;
            Up = Up & (ind(d,:)+abs(Q(d,q))<L(d)) & (ind(d,:)-abs(Q(d,q))>-1) ;
            
    end
end
selectInd = [Dwn1 ; Dwn2 ; Up] ;
selectInd = reshape(selectInd',[L 3]) ; 
selectInd = reshape(selectInd,[],3) ;

xx = ii(:,:,1) ;
yy = ii(:,:,2) ;
    for s = 1:3
        ax(s) = mysubplot(1,3,s) ;
            plot(xx(:),yy(:),'.','markersize',5,'linewidth',1,'color',[1 1 1]*0.)
            plot(xx(selectInd(:,s)),yy(selectInd(:,s)),marker,'color',colors(q,:),'markersize',markersize,'linewidth',linewidth) ;
            ax(s).XTick = [] ; ax(s).YTick = [] ;
            axis equal
            axis tight
            axis off
            ax(s).OuterPosition = ax(s).OuterPosition - 0.07*[0 -0.85 0 1.85] ;
            %xlabel(['(',char(96+q),')', ' $\vect q = [',num2str(Q(1,q)),'\;\;',num2str(Q(2,q)),']$' ],'visible','on') ;
            %ax.OuterPosition = [(q-1)*1/nQ 0 1/nQ 1] ;
    end
    % Make the second axis closer to the first
        ax(2).OuterPosition(1) = ax(2).OuterPosition(1)-.04 ;
    % Annotations axes
        ax(4) = axes('position',[0 0 1 1]) ;
            axis off
            text(sum([ax(1).Position(1),ax(2).Position([1,3])])/2,.5,...
                '\large{$+$}'...
                ...,'color',colors(q,:)...
                )
            text(sum([ax(2).Position(1),ax(3).Position([1,3])])/2,.5,...
                ...['\Large{$ = \mat \Pi^{\begin{bmatrix}',num2str(Q(1,q)),'\\',num2str(Q(2,q)),'\end{bmatrix}}$}']...
                '\large{$ = 2 \hsp \mat \Pi^{\vect q} \hsp \times$}'...
                ...,'color',colors(q,:)...
                )

% SAVE THE TIKZ FIGURE 
    cleanfigure ;
    figPos = get(gcf,'position') ;
    normHeight = figPos(4)/figPos(3) ;
    sufix = ['_',num2str(q)] ;
    %sufix = ['fun_',fun{q}(1,:),'_',fun{q}(2,:),'_Q_',num2str(Q(1,q)),'_',num2str(Q(2,q))] ;
    matlab2tikz(['InvRotGrids',sufix,'.tikz']...
                ,'width','\figwidth'...
                ,'height',[num2str(normHeight,4),'\figwidth']...
                ...,'noSize',true...
                ...,'parseStringsAsMath',true...
                ) ; 
end   

% Create the overall file
    Str = [] ;
    for q=1:nQ
        Str = [Str,'\begin{subfigure}{\textwidth}',char(10)] ;
        Str = [Str,'\centering',char(10)] ;
        Str = [Str,char(10)] ;
        Str = [Str,char(10)] ;
        
        fName = ['InvRotGrids','_',num2str(q),'.tikz'] ;
        fID = fopen(fName,'r') ;
        StrFile = fread(fID) ;
        fclose(fID) ;
        delete(fName) ;
        Str = [Str,char(StrFile(:)'),char(10)] ;
        
        Str = [Str,char(10)] ;
        Str = [Str,char(10)] ;
        
        de = ismember(num2cell(fun{q},2),{'exp'})' ;
        dc = ismember(num2cell(fun{q},2),{'cos'})' ;
        if sum(de)==0 
            strDe = '\emptyset';
        else
            dimsDe = double(de).*[1:D] ;
            strDe = regexprep(mat2str(dimsDe(dimsDe~=0)),' ',' \\;\\; ') ;
        end
        if sum(dc)==0 
            strDc = '\emptyset';
        else
            dimsDc = double(dc).*[1:D] ;
            strDc = regexprep(mat2str(dimsDc(dimsDc~=0)),' ',' \\;\\; ') ;
        end
        
        %Str = [Str,'\caption{$\vect q = [',num2str(Q(1,q)),' \;\; ',num2str(Q(2,q)),']$, $(D\idr{e},D\idr{c}) = (',num2str(De),',',num2str(Dc),')$}',char(10)] ;
        Str = [Str,'\caption{$\vect q = [',num2str(Q(1,q)),' \;\; ',num2str(Q(2,q)),']$, $\vect d\idr{e} = ',strDe,'$, $\vect d\idr{c} = ',strDc,'$.}',char(10)] ;
        Str = [Str,'\end{subfigure}',char(10)] ;
        Str = [Str,'\hfill \vspace{\Vspacing}',char(10)] ;
        
        Str = [Str,char(10)] ;
        Str = [Str,char(10)] ;
        Str = [Str,char(10)] ;
    end
    fID = fopen('InvRotGridsSamples.tex','w') ;
    fwrite(fID,Str) ;
    fclose(fID) ;

close all ;
