% Read header of a VAC data file

fileerror=0;
if asciifile(ifile)
   headline=trim(fgetl(fid));
   if ~isstr(headline);fileerror=1;return;end;
   tmp=sscanf(fgetl(fid),'%d %f %d %d %d',5)';
   it=tmp(1); time=tmp(2); ndim=tmp(3); neqpar=tmp(4); nw=tmp(5); clear tmp;
   gencoord=ndim<0; ndim=abs(ndim);
   nx=sscanf(fgetl(fid),'%d',ndim)';
   eqpar=sscanf(fgetl(fid),'%f',neqpar)';
   Variables=trim(fgetl(fid));
else
   [tmp,ntmp]=fread(fid,4);
   if ntmp<4;fileerror=1;return;end;
   headline=trim(setstr(fread(fid,79)'));
   fread(fid,4);

   fread(fid,4);
   it=fread(fid,1,'int32'); time=fread(fid,1,'float64');
   ndim=fread(fid,1,'int32');
   gencoord=ndim<0; ndim=abs(ndim);
   neqpar=fread(fid,1,'int32'); nw=fread(fid,1,'int32');
   fread(fid,4);

   fread(fid,4);
   nx=fread(fid,ndim,'int32');
   fread(fid,4);

   fread(fid,4);
   eqpar=fread(fid,neqpar,'float64');
   fread(fid,4);

   fread(fid,4);
   Variables=trim(setstr(fread(fid,79)'));
   fread(fid,4);
end

% Extract physics from headline if it is defined
i=length(headline);
while i>1 & headline(i)~='_'; i=i-1; end;
if headline(i)=='_'
   physics=headline(i+1:length(headline));
   headline=headline(1:i-1);
end
% Extraxt number of vector components NDIR from last character of physics
% and extract phys without the number of dimesions and components
if ~isempty(physics)
   ndir=str2num(physics(length(physics)));
   phys=physics(1:length(physics)-2);
end

% Extract info from names of Variables
[variables,ntmp]=str2arr(Variables);
xnames=variables(1:ndim,:);
wnames=variables(ndim+1:ndim+nw,:);

% It can optionally contain the names of the equation parameters
if ntmp==ndim+nw+neqpar
   eqparnames=variables(ndim+nw+1:ntmp,:);
   for ieqpar=1:neqpar
      eval([eqparnames(ieqpar,:) '= eqpar(ieqpar);']);
   end
end
clear ntmp;

nxs=prod(nx);
if ndim==1
   nx1=nx;
   nx2=1;
   nx3=1;
end
if ndim==2
   nx1=nx(1);
   nx2=nx(2);
   nx3=1;
end
if ndim==3
   nx1=nx(1);
   nx2=nx(2);
   nx3=nx(3);
end
