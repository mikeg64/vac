function [xreg_out,yreg_out]=triangulate(xirr_in,yirr_in,nx,ny)
% Set the global NEI1,NEI2,NEI3,DST1,DST2,DST3 variables for interpolate
% the xirr_in, yirr_in vectors contain the coordinates of the irregular points,
% while nx, ny gives the size of the regular grid onto which we interpolate.
% The output vectors xreg_out,yreg_out contain the regular coordinates.

% Neighbour indices and distances
global NEI1 NEI2 NEI3 DST1 DST2 DST3 DSTMAX XREG YREG XIRR YIRR;

disp('Triangulate...');

XIRR=xirr_in; YIRR=yirr_in;

% Set up the regular grid (normalized would be better!!!)
xmax=max(max(XIRR)); ymax=max(max(YIRR));
xmin=min(min(XIRR)); ymin=min(min(YIRR));
dx=(xmax-xmin)/(nx-1);   dy=(ymax-ymin)/(ny-1);
[YREG,XREG]=meshgrid(ymin:dy:ymax,xmin:dx:xmax);

DSTMAX=2*(xmax-xmin+ymax-ymin)+1;

NEI1=zeros(nx,ny);NEI2=zeros(nx,ny);NEI3=zeros(nx,ny);
DST1=DSTMAX*ones(nx,ny); DST2=DST1; DST3=DST1;

% Indices of closest regular points left and below the irregular ones
indirr=fix((XIRR-xmin)/dx-0.00001)+nx*fix((YIRR-ymin)/dy-0.00001)+1;

% Put neighbour info into regular points surrounding irregular points
% The double loop is for the four corners of the enclosing regular cell
for dix=0:1
   for diy=0:1
      ind0=indirr+dix+nx*diy; 
      dd=abs(XIRR-XREG(ind0))+abs(YIRR-YREG(ind0));
      % We allow four irregular points in a regular cell, 
      % each gets a chance to propagate into a corner
      for i=1:4
         iind=find(DST1(ind0)>dd);
         if ~isempty(iind)
            ind=ind0(iind);
            NEI1(ind)=iind;
            DST1(ind)=dd(iind);
         end
      end
   end
end

% Index arrays for fast shifting in the four directions
[Aiy,Aix]=meshgrid(1:ny,1:nx);
jxiy=Aix+(Aix<nx)+nx*(Aiy-1);hxiy=Aix-(Aix>1)+nx*(Aiy-1);
ixjy=Aix+nx*(Aiy-1+(Aiy<ny));ixhy=Aix+nx*(Aiy-1-(Aiy>1));

% Iteration to fill up NEI1,NEI2,...,DST3 arrays
% We are done if NEI3 is filled, we fail if 100 iterations are not enough
count=100;
while count>0 & ~all(all(NEI3))
   count=count-1;
   ind0=find(NEI3==0 & NEI2>0);
   if ~isempty(ind0)
      trinei(3,ind0,jxiy);
      trinei(3,ind0,hxiy);
      trinei(3,ind0,ixjy);
      trinei(3,ind0,ixhy);
   end
   ind0=find(NEI2==0 & NEI1>0);
   if ~isempty(ind0)
      trinei(2,ind0,jxiy);
      trinei(2,ind0,hxiy);
      trinei(2,ind0,ixjy);
      trinei(2,ind0,ixhy);
   end
   ind0=find(NEI1==0);
   if ~isempty(ind0)
      trinei(1,ind0,jxiy);
      trinei(1,ind0,hxiy);
      trinei(1,ind0,ixjy);
      trinei(1,ind0,ixhy);
   end
end
if ~count
   disp('Did not find neighbors in 100 iterations!');
end

% Transform distances to weights and put result back into DST1,DST2,DST3
dst1=DST2.*DST3; dst2=DST1.*DST3; dst3=DST1.*DST2; dd=dst1+dst2+dst3;
DST1=dst1./dd; DST2=dst2./dd; DST3=dst3./dd;

% Output regular grid
xreg_out=XREG; yreg_out=YREG;

% Clear unused globals
clear global DSTMAX XREG YREG XIRR YIRR;
