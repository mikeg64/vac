function freg=interpolate(f)
% interpolate matrix f (on irregular grid) into freg (on regular grid)
% using NEI1,NEI2,NEI3 neighbour indices and DST1,DST2,DST3 weights
% The indices and weights should be calculated by regular_grid

global NEI1 NEI2 NEI3 DST1 DST2 DST3

if isempty(NEI1)
   disp('Call triangulate first!')
   return
end

freg=DST1.*f(NEI1)+DST2.*f(NEI2)+DST3.*f(NEI3);

