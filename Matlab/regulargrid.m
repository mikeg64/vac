% Produce regulargrid and interpolte variables onto it
% The size of the regular grid is nxreg, the variables are taken from wnames

if isempty(nxreg)
   Transform='regular';
   read_transform_par;
end

% Check whether nxreg or X has changed since last triangulate
newX=1;
if ~isempty(nxregold)
   if all(nxregold==nxreg)
      if all(size(Xold)==size(X))
         newX= any(any(Xold~=X));
      end
   end
end

if newX
   Xold=X;
   nxregold=nxreg;
   [xreg,yreg]=triangulate(xx,yy,nxreg(1),nxreg(2));
end

% Redefine coordinates and interpolate variables
xx=xreg;
yy=yreg;
eval([xnames(1,:) '=xreg;']);
eval([xnames(2,:) '=yreg;']);
for iw=1:nw
   tmp=wnames(iw,:);
   eval([tmp '=interpolate(' tmp ');']);
end
