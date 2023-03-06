% Plot functions defined by funcs according to plotmodes

if plotdim>2
   disp('plot_func.m cannot plot 3D yet');
   return
end

if isempty(cut)
   x_=xx;
   if ndim==2;y_=yy;end;
else
   x_=eval(['xx(' cut ')']);
   if ndim==2;y_=eval(['yy(' cut ')']);end;
end
if size(x_,2)==1
   xlimits=[min(x_) max(x_)];
elseif size(x_,1)==1
   xlimits=[min(y_) max(y_)];
   x_=y_;
else
   xlimits=[min2(x_) max2(x_) min2(y_) max2(y_)];
end

for ifunc=1:nfunc
   if multidir
      plotix=floor(multi0/multiy)+1;
      plotiy=multi0-multiy*floor(multi0/multiy)+1;
      subplot(multiy,multix,(plotiy-1)*multix+plotix);
   else
      plotix=multi0-multix*floor(multi0/multix)+1;
      plotiy=floor(multi0/multix)+1;
      subplot(multiy,multix,multi0+1);
   end;
   multi0=multi0+1;

   fun=funcs(ifunc,:);
   get_func;
   plotmod=trim(plotmodes(ifunc,:));
   flimits=[fmin(ifunc) fmax(ifunc)];
   if flimits(1)==flimits(2)
      flimits(1)=flimits(1)-1; flimits(2)=flimits(2)+1;
   end;
   limits=[xlimits flimits]; 

   if plotdim==1
      if isempty(findstr(['/' plotmod '/'],'/plot/loglog/semilogx/semilogy/'))
         if isfirstpict
            disp(['Warning: using default plotmode plot instead of ''' ...
                   plotmod ''' ']);
         end
         plotmod='plot';
      end
      eval([plotmod '(x_,f1)']);
      axis(limits);
   end
   if plotdim==2
   if strcmp(plotmod,'surf')
      surf(x_,y_,f1);
      view(View);
      axis(limits);
   elseif strcmp(plotmod,'surfc')
      surf(x_(:,1),y_(1,:),f1);
      view(View);
      axis(limits);
   elseif strcmp(plotmod,'surfl')
      surfl(x_,y_,f1);
      view(View);
      axis(limits);
   elseif strcmp(plotmod,'pcolor')
      pcolor(x_,y_,f1);
      axis(xlimits);
   elseif strcmp(plotmod,'contour')
      contour(x_(:,1),y_(1,:),f1',Contourlevel,Contourstyle);
      axis(xlimits);
   elseif strcmp(plotmod,'quiver')
      if NF==2
         quiver(x_,y_,f1,f2,Quiverscale);
         axis(xlimits);
      elseif isfirstpict
         disp('Error: quiver requires 2 function components!');
      end
   elseif strcmp(plotmod,'pcolor_contour') | strcmp(plotmod,'pcolorc')
      hold on;
      pcolor(x_,y_,f1);
      if NF==1;f2=f1;end;
      contour(x_(:,1),y_(1,:),f2',Contourlevel,Contourstyle);
      axis(xlimits);
   elseif strcmp(plotmod,'contour_quiver')
      hold on;
      contour(x_(:,1),y_(1,:),f1',Contourlevel,Contourstyle);
      if NF==3
         quiver(x_,y_,f2,f3,Quiverscale);
      elseif isfirstpict
         disp('Error: contour_quiver requires 3 function components!');
      end
      axis(xlimits);
   elseif strcmp(plotmod,'pcolor_quiver')
      hold on;
      pcolor(x_,y_,f1);
      if NF==3
         quiver(x_,y_,f2,f3,Quiverscale);
      elseif isfirstpict
         disp('Error: pcolor_quiver requires 3 function components!');
      end
      axis(xlimits);
   else
      if ~strcmp(plotmod,'mesh') & isfirstpict
         disp(['Warning: using default plotmode mesh instead of ''' ...
                plotmod ''' ']);
      end
      mesh(x_,y_,f1);
      view(View);
      axis(limits);
   end
   end % of plotdim==2

   if size(x_,2)==1
      xlabel(xnames(1));
   elseif size(x_,1)==1
      xlabel(xnames(2));
   else
      xlabel(xnames(1));
      ylabel(xnames(2));
      eval(['shading ' Shadings(ifunc,:)]);
      if ~strcmp(plotmode,'surfl')
         caxis(flimits);
         if Colorbar(ifunc);colorbar;end;
      end
   end

   tmp=trim(plottitles(ifunc,:));
   if strcmp(tmp,'default')
      title(funcs(ifunc,:));
   else
      title(tmp);
   end;
   hold off;
end

if Bottomline>0 | Headerline>0
   if npict1==1
      axes('position',[(ifile-1)/nfile 0 1/nfile 1], ...
           'visible','off','nextplot','new');
      % The nextplot new may not be needed at all or maybe HandleVisibility
      % axes('position',[(ifile-1)/nfile 0 1/nfile 1],'visible','off');
      if Headerline>0
         tmp=headline;
         if Headerline>1;tmp=[tmp ' (nx=' sprintf('%d ',nx) ')']; end;
         text(0.01,0.99,tmp);
      end
   else
      if Headerline>0 & ipict1==0
         axes('position',[(ifile-1)/multix (multiy-1)/multiy ...
              1/multix 1/multiy],'visible','off','nextplot','new');
         % The nextplot new may not be needed at all or maybe HandleVisibility 
         % axes('position',[(ifile-1)/multix (multiy-1)/multiy ...
         %       1/multix 1/multiy],'visible','off','HandleVisibility','on');
         tmp=headline;
         if Headerline>1;tmp=[tmp ' (nx=' sprintf('%d ',nx) ')']; end;
         text(0.01,0.99,tmp);
      end;
      if Bottomline>0
         % The nextplot new may not be needed at all or maybe HandleVisibility
         axes('position',[(plotix-1)/multix (multiy-plotiy)/multiy...
              1/multix 1/multiy],'visible','off','nextplot','new');
      end;
   end;

   if Bottomline>0
      tmp=sprintf('time=%9.3g',time);
      if Bottomline>1;tmp=[tmp sprintf(' it=%5d',it)]; end;
      if Bottomline>2;tmp=[tmp ' nx=' sprintf('%d ',nx)]; end;
      text(0.01,0.02,tmp);
   end;
end;
