function [dst,nei]=triset(dst_in,nei_in,dd,nn,ind0,iind);
% Set distance dst and neighbour index nei of ind0[iind] from nn and dd.

dst=dst_in;
nei=nei_in;
if ~isempty(iind)
   ind=ind0(iind);
   dst(ind)=dd(iind);
   nei(ind)=nn(iind);
end
