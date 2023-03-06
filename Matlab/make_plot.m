for j=1:nplots

   subplot(subpl(1),subpl(2),j)
   
   if(findstr(plotmodemat(j,:),'mesh')~='')
      funcact=funcmat(j,:);
      get_plot
      mesh(x(ii:ie,ji:je),y(ii:ie,ji:je),plv(ii:ie,ji:je))
      view(az,el)
      axis(ranges(j,:))
      xlabel('x')
      ylabel('y')
      if(colorbarmesh=='y')
         colorbar
      end
   end

   if(findstr(plotmodemat(j,:),'mereg')~='')
      funcact=funcmat(j,:);
      get_plot
      interpolfunc
      mesh(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)')
      view(az,el)
      axis([xi(ii),xi(ie),yi(ji),yi(je),ranges(j,5),ranges(j,6)])
      xlabel('x')
      ylabel('y')
      if(colorbarmesh=='y')
         colorbar
      end
   end

   if(findstr(plotmodemat(j,:),'quivcona')~='')
      funcact=funcmat(j,7:size(funcmat,2));
      get_plot
      interpolfunc
      contour(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)',ncontour)
      xlabel('x')
      ylabel('y')
      hold
      funcact=funcmat(j,1:2);
      get_plot
      interpolfuncnolim
      plv1=plv;
      funcact=funcmat(j,4:5);
      get_plot
      interpolfuncnolim
      plv2=plv;
      quiver(xi(ii:ie),yi(ji:je),plv1(ii:ie,ji:je)',plv2(ii:ie,ji:je)',quivers)
      if(colorbarquivcont=='y')
         colorbar
      end
      hold off
   end

   if(findstr(plotmodemat(j,:),'quivcont')~='')
      funcact=funcmat(j,7:size(funcmat,2));
      get_plot
      interpolfunc
      vcont=linspace(ranges(j,5),ranges(j,6),ncontour);
      contour(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)',vcont)
      xlabel('x')
      ylabel('y')
      hold
      funcact=funcmat(j,1:2);
      get_plot
      interpolfuncnolim
      plv1=plv;
      funcact=funcmat(j,4:5);
      get_plot
      interpolfuncnolim
      plv2=plv;
      quiver(xi(ii:ie),yi(ji:je),plv1(ii:ie,ji:je)',plv2(ii:ie,ji:je)',quivers)
      if(colorbarquivcont=='y')
         colorbar
      end
      hold off
   end
   
   if(findstr(plotmodemat(j,:),'pcolcona')~='')
      funcact=funcmat(j,1:6);
      get_plot
      interpolfuncnolim
      colormap('cool')
      pcolor(x(ii:ie,ji:je),y(ii:ie,ji:je),plv(ii:ie,ji:je))
      if(pcolorflat=='y')
         shading flat
      end
      axis(ranges(j,1:4))
      xlabel('x')
      ylabel('y')
      hold on
      funcact=funcmat(j,7:size(funcmat,2));
      get_plot
      interpolfunc
      colormap('hsv')
      contour(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)',ncontour)
      xlabel('x')
      ylabel('y')
      hold off
   end

   if(findstr(plotmodemat(j,:),'quivpcol')~='')
      funcact=funcmat(j,7:size(funcmat,2));
      get_plot
      pcolor(x(ii:ie,ji:je),y(ii:ie,ji:je),plv(ii:ie,ji:je))
      if(pcolorflat=='y')
         shading flat
      end
      axis(ranges(j,1:4))
      xlabel('x')
      ylabel('y')
      hold
      funcact=funcmat(j,1:2);
      get_plot
      plv1=plv;
      funcact=funcmat(j,4:5);
      get_plot
      plv2=plv;
      quiver(x(ii:ie,ji:je),y(ii:ie,ji:je),plv1(ii:ie,ji:je),plv2(ii:ie,ji:je),quivers)
      if(colorbarquivpcol=='y')
         colorbar
      end
      hold off
   end

   if(findstr(plotmodemat(j,:),'pcolor')~='')
      funcact=funcmat(j,:);
      get_plot
      pcolor(x(ii:ie,ji:je),y(ii:ie,ji:je),plv(ii:ie,ji:je))
      if(pcolorflat=='y')
         shading flat
      end
      axis(ranges(j,1:4))
      xlabel('x')
      ylabel('y')
      if(colorbarpcolor=='y')
         colorbar
      end
   end
   
   if(findstr(plotmodemat(j,:),'contour')~='')
      funcact=funcmat(j,:);
      get_plot
      interpolfunc
      vcont=linspace(ranges(j,5),ranges(j,6),ncontour);
      contour(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)',vcont)
      xlabel('x')
      ylabel('y')
      if(colorbarcontour=='y')
         colorbar
      end
   end

   if(findstr(plotmodemat(j,:),'conta')~='')
      funcact=funcmat(j,:);
      get_plot
      interpolfunc
      contour(xi(ii:ie),yi(ji:je),plv(ii:ie,ji:je)',ncontour)
      xlabel('x')
      ylabel('y')
      if(colorbarcontour=='y')
         colorbar
      end

   end
   
   if(findstr(plotmodemat(j,:),'quiver')~='')
      funcact=funcmat(j,1:2);
      get_plot
      plv1=plv;
      funcact=funcmat(j,4:5);
      get_plot
      plv2=plv;
      quiver(x(ii:ie,ji:je),y(ii:ie,ji:je),plv1(ii:ie,ji:je),plv2(ii:ie,ji:je),quivers)
      axis(ranges(j,1:4))
      xlabel('x')
      ylabel('y')
   end

   if(j==1)
      title(['time=' num2str(t) ', iter=' num2str(step) ', ' date  ', ' dirname ' / ' funcmat(j,:)])
   else
      title(funcmat(j,:))
   end

end
