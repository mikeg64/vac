function trinei(inei,ind0,shift)
% Propagate inei-th neighbour of ind0 from direction given by shift

global NEI1 NEI2 NEI3 DST1 DST2 DST3 DSTMAX XREG YREG XIRR YIRR;

jnd0=shift(ind0);

% Starting with the inei-th neighbour we go down to the 1st neighbours
% of jnd0 to find the new inei-th neighbour of ind0.
for i=inei:-1:1
   if i==1
      nn=NEI1(jnd0);
   elseif i==2
      nn=NEI2(jnd0);
   else
      nn=NEI3(jnd0);
   end
   % Get indices and values of neighbours which are already known (>0)
   [inn1,dummy,nn1]=find(nn);
   % Set distances to DSTMAX for all then set distances for these
   dd=DSTMAX+0*ind0;
   if ~isempty(inn1)
     dd(inn1)=abs(XREG(ind0(inn1))-XIRR(nn1))+abs(YREG(ind0(inn1))-YIRR(nn1));
   end
   % By comparing the dd distances with the current distances 
   % propagate neighbour and distance info from nn and dd
   if inei==1
      % Propagate 1st neighbour of ind0 from nn and dd
      iind=find(DST1(ind0)>dd);
      [DST1,NEI1]=triset(DST1,NEI1,dd,nn,ind0,iind);
   elseif inei==2
      % Propagate 2nd neighbour of ind0 from nn and dd
      iind=find(DST2(ind0)>dd & NEI1(ind0)~=nn);
      [DST2,NEI2]=triset(DST2,NEI2,dd,nn,ind0,iind);
   else
      % Propagate 3rd neighbour of ind0 from nn and dd
      iind=find(DST3(ind0)>dd & NEI1(ind0)~=nn & NEI2(ind0)~=nn);
      [DST3,NEI3]=triset(DST3,NEI3,dd,nn,ind0,iind);
   end
end