% Read a snapshot from a VAC data file (binary or ascii).
% First X and w are read, (X is capitalized to avoid possible conflict with 
% the name of its first component). Transform X and w according to Transform
% if generalized coordinates are found.
% Extract the variables given in the xwname string array read by get_head.

X=zeros(nxs,ndim);
w=zeros(nxs,nw);
if asciifile(ifile)
   xw=fscanf(fid,'%f',[ndim+nw,nxs]);
   X=xw(1:ndim,:)';
   w=xw(ndim+1:nw+ndim,:)';
   clear xw;
   fgetl(fid);
else
   fread(fid,4);
   for idim=1:ndim
      X(:,idim)=fread(fid,nxs,'float64');
   end
   fread(fid,4);
   for iw=1:nw
      fread(fid,4);
      w(:,iw)=fread(fid,nxs,'float64');
      fread(fid,4);
   end
end

% extract variables from X into variables named after the strings in xnames
for idim=1:ndim
  if ndim==2
     tmp=reshape(X(:,idim),nx1,nx2);
  else
     tmp=X(:,idim);
  end
  eval([xnames(idim,:),'=tmp;']);
end
if ndim==1
   xx=X;
elseif ndim==2
   xx=reshape(X(:,1),nx1,nx2);
   yy=reshape(X(:,2),nx1,nx2);
end

% extract variables from w into variables named after the strings in wnames
for iw=1:nw
  if ndim==2
     tmp=reshape(w(:,iw),nx1,nx2);
  else
     tmp=w(:,iw);
  end
  eval([wnames(iw,:),'=tmp;']);
end
clear tmp;

if gencoord & ndim==2 
   if strcmp(Transform,'polar')
      polargrid;
   elseif strcmp(Transform,'regular')
      regulargrid;
   end
end
