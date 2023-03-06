% Read the plotting parameters for plotfunc or animate

func     =askstr('func(s) (rho m1;m2 ...)       ',func);
[funcs,nfunc]=str2arr(func);
if ~isempty(cut); disp(['cut = ' cut]); end;

autorange=asknum('autorange(s) (1/0 or [1 0 1]) ',autorange,nfunc);

if isempty(cut)
   plotdim=ndim;
elseif ndim==1
   plotdim=1;
elseif ndim==2
   if gencoord & strcmp(Transform,'regular')
      tmp=zeros(nxreg(1),nxreg(2));
   else
      tmp=zeros(nx1,nx2);
   end
   tmp=eval(['tmp(' cut ')']);
   if size(tmp,1)==1 | size(tmp,2)==1
      plotdim=1;
   else
      plotdim=2;
   end
   clear tmp;
else
   display('plotdim=2 is a guess for ndim>2');
   plotdim=2
end
if plotdim==1
   disp('1D plotmodes: plot/loglog/semilogx/semilogy');
elseif plotdim==2
   disp(['2D plotmodes for 1 variable'...
        ' : mesh/surf/surfc/surfl/contour/pcolor/pcolorc']);
   disp( '             for 2 variables: quiver/pcolor_contour');
   disp( '             for 3 variables: contour_quiver/pcolor_quiver');
else
   disp('No 3D plotmodes yet...');
   return
end

plotmode  =askstr('plotmode(s)                   ',plotmode);
plotmodes =str2arr(plotmode,nfunc);
plottitle =askstr('plottitle(s) (e.g. ''B [G];J'') ',plottitle);
plottitles=comma2arr(plottitle,nfunc);

if plotdim>1
   Shading=  askstr('Shading(s) flat/faceted/interp',Shading);
   Shadings= str2arr(Shading,nfunc);
   Colorbar =asknum('Colorbar(s) (1/0 or [1 0 1])  ',Colorbar,nfunc);
   fprintf(['View=[%g %g] '],View);
 fprintf(['Contourlevel=%d Contourstyle=''' Contourstyle ...
          ''' Quiverscale=%g\n'],Contourlevel,Quiverscale);
end

fprintnum('multiplot=',multiplot);
fprintf([' Bottomline=%d Headerline=%d\n'],Bottomline,Headerline);
