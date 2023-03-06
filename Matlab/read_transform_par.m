% Read the transformation parameters for irregular grids

Transform=askstr('Transform (regular/polar/none)',Transform);
if strcmp(Transform,'regular')
   nxreg=asknum('nxreg (eg. 100 or [50 60])',nxreg,2);
elseif strcmp(Transform,'polar')
   disp('Radial and phi component(s) of vector variables among');
   fprintf(          'wnames                           = ');
   for iw=1:nw; fprintf([trim(wnames(iw,:)),' ']); end; 
   if ~isempty(polar_r); fprintf('\n'); end;
   polar_r   =askstr('polar_r(s)   (eg. mr br     ...)',polar_r);
   polar_phi =askstr('polar_phi(s) (eg. mphi bphi ...)',polar_phi);
   polar_rs  =str2arr(polar_r);
   polar_phis=str2arr(polar_phi);
elseif ~strcmp(Transform,'none')
   disp(['Incorrect value ' Transform ' is replaced by ''none''']);
   Transform='none'
end
