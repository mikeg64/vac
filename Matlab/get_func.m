% Calculate f1 (or f1,f2..) from the function name string fun.
% fun may contain more variables separated by semicolons, e.g.
%
% m1./rho;m2./rho
%
% In this case f1=m1./rho, f2=m2./rho is returned.
%
% For each variable we first check if it is defined as a function 
% (e.g. v1 or v2), otherwise the string is assumed to be an expression
% of the conserved variables defined by get_body and it is evaluated with EVAL.
%
% If the "cut" string is defined, the variables are cut, e.g. cut=':,10'
% will take a cut along the first dimension at the 10-th row.

[funs,NF]=comma2arr(fun);

for i=1:NF
   fun=trim(funs(i,:));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Put new function definitions here below!                           %
   % Use "fun", "physics", "phys", and "ndir" as switches.              %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if strcmp(fun,'v1')
      tmp=m1./rho;
   elseif strcmp(fun,'v2')
      tmp=m2./rho;
   elseif strcmp(fun,'v3')
      tmp=m3./rho;
   elseif strcmp(fun,'pth') | strcmp(fun,'T') | strcmp(fun,'csound') | ...
          strcmp(fun,'mach')
      if ~(strcmp(phys,'hdadiab') | strcmp(phys,'hd')| ...
           strcmp(phys,'hd')      | strcmp(phys,'mhd'))
         disp(['physics=''' physics ''' is not known for function ''' ...
                fun ''' !']);
         tmp=0*xx;
      else
         % Ram pressure=rho*v^2 for mach number or HD and MHD
         tmp_ram=m1.^2;
         if ndir>1; tmp_ram=tmp_ram+m2.^2; end;
         if ndir>2; tmp_ram=tmp_ram+m3.^2; end;
         tmp_ram=tmp_ram./rho;

         % Put thermal pressure into tmp
         if strcmp(phys,'hdadiab')
            tmp=eqpar(2)*rho.^eqpar(1);
         elseif strcmp(phys,'mhdiso')
            tmp=eqpar(2)*rho;
         else
            % Add twice the magnetic energy b^2
            if strcmp(phys,'mhd')
               tmp=tmp_ram+b1.^2;
               if ndir>1; tmp=tmp+b2.^2; end;
               if ndir>2; tmp=tmp+b3.^2; end;
            end
            tmp=(eqpar(1)-1)*(e-tmp/2);
         end;
         % Temperature=pth/rho
         % Sound speed=sqrt(gamma*pth/rho)
         % Mach number=v/csound=sqrt((rho*v^2)/(gamma*pth))
         if strcmp(fun,'T');          tmp=tmp./rho;
         elseif strcmp(fun,'csound'); tmp=sqrt(eqpar(1)*tmp./rho);
         elseif strcmp(fun,'mach');   tmp=sqrt(tmp_ram./(eqpar(1)*tmp));
         end;
         clear tmp_ram;
      end;
   elseif strcmp(fun,'curlv')
      tmp=gradx(m2./rho,xx)-grady(m1./rho,yy);
   elseif strcmp(fun,'divb')
      tmp=gradx(b1,xx)+grady(b2,yy);
   elseif strcmp(fun,'divb4')
      tmp=grad4x(b1,xx)+grad4y(b2,yy);
   elseif strcmp(fun,'divb_rz')
      tmp=gradx(xx.*b1,xx)./xx+grady(b2,yy);
   elseif strcmp(fun,'divb_rp')
      tmp=gradx(xx.*b1,xx)./xx+grady(b2,yy)./xx;
   elseif strcmp(fun,'j')
      tmp=gradx(b2,xx)-grady(b1,yy);
   elseif strcmp(fun,'j_rz')
      tmp=grady(b1,yy)-gradx(b2,xx);
   elseif strcmp(fun,'j_rp')
      tmp=gradx(b2,xx)-grady(b1,yy)./xx;
   elseif strcmp(fun,'A')+strcmp(fun,'A_r')+strcmp(fun,'AA')+strcmp(fun,'AA_r')
      % Calculate vector potential from b1=dA/dy b2=-dA/dx
      bx=b1;
      by=b2;
      if strcmp(fun,'A_r')+strcmp(fun,'AA_r')
         % For axial symmetry Cartesian_curl(r*A)=(B_r*r,B_z*r)
         bx=xx*bx;
         by=xx*by;
      end
      tmp=zeros(nx1,nx2);
      if strcmp(fun,'A')+strcmp(fun,'A_r')
         % Integrate along the first row
         for i1=2:nx1; 
            tmp(i1,1)=tmp(i1-1,1) ...
               +(by(i1,1)+by(i1-1,1))*(xx(i1,1)-xx(i1-1,1))*0.5 ...
               -(bx(i1,1)+bx(i1-1,1))*(yy(i1,1)-yy(i1-1,1))*0.5; 
         end
         % Integrate all columns vertically
         for i2=2:nx2;
            tmp(:,i2)=tmp(:,i2-1) ...
               +(by(:,i2)+by(:,i2-1)).*(xx(:,i2)-xx(:,i2-1))*0.5 ...
               -(bx(:,i2)+bx(:,i2-1)).*(yy(:,i2)-yy(:,i2-1))*0.5;
         end
      else
         % Integrate first column vertically
         for i2=2:nx2
            tmp(1,i2)=tmp(1,i2-1) ...
                -(bx(1,i2)+bx(1,i2-1))*(yy(1,i2)-yy(1,i2-1))*0.5 ...
                +(by(1,i2)+by(1,i2-1))*(xx(1,i2)-xx(1,i2-1))*0.5;
         end
         % Integrate all rows
         for i1=2:nx1 
           tmp(i1,:)=tmp(i1-1,:) ...
                -(bx(i1,:)+bx(i1-1,:)).*(yy(i1,:)-yy(i1-1,:))*0.5 ...
                +(by(i1,:)+by(i1-1,:)).*(xx(i1,:)-xx(i1-1,:))*0.5;
         end
      end
   else
      % The default assumption is an expression of the conserved variables.
      tmp=eval(funs(i,:));
   end
   % Cut out the part to be plotted and put the result in f1 (f2,...)
   if isempty(cut)
      eval(['f' num2str(i) '=tmp;']);
   else
      eval(['f' num2str(i) '=tmp(' cut ');']);
   end;
end;
