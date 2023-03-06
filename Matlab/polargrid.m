% Transform the polar coordinates and the polar vector variables
% The coordinate names are stored in xnames. The generic xx and yy are
% also transformed.
% The vector variable names are stored in polar_rs and polar_phis.

if ndim~=2
   disp('polargrid works for 2D arrays only!');
   return;
end
if isempty(polar_r) | isempty(polar_phi)
   Transform='polar'; read_transform_par;
end;

% Calculate R and PHI from the original X and Y read from the data file
tmp1=sqrt(xx.^2+yy.^2);
tmp2=atan2(yy,xx);
xx=tmp1;
yy=tmp2;
% Remove 2*pi jumps from PHI=yy
for i=2:nx2;if yy(1,i-1)>yy(1,i);yy(:,i)=yy(:,i)+2*pi;end;end;
% Put R and PHI into the coordinate variables defined by xnames
eval([xnames(1,:) '=xx;']);
eval([xnames(2,:) '=yy;']);

% Rotate polar variables
for i=1:size(polar_r,1)
   tmp1= eval(polar_rs(i,:)).*cos(yy)+eval(polar_phis(i,:)).*sin(yy);
   tmp2=-eval(polar_rs(i,:)).*sin(yy)+eval(polar_phis(i,:)).*cos(yy);
   eval([polar_rs(i,:)   '=tmp1;']);
   eval([polar_phis(i,:) '=tmp2;']);
end

clear tmp1 tmp2;
