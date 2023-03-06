function [isascii,filesize,pictsize,fileerror]=studyfile(filenames)
% Determine whether the files are ascii (1) or binary (0) based on 
% the first 4 characters. Determine the filesize and the size of one snapshot.
% Fileerror is set to 1 if any of the files could not be opened. 
% Usage:
%    [isascii,filesize,pictsize,fileerror]=studyfile(filenames)

nfile=size(filenames,1);
isascii=zeros(1,nfile);
filesize=zeros(1,nfile);
pictsize=zeros(1,nfile);
fileerror=0;

for ifile=1:nfile
   fid=fopen(trim(filenames(ifile,:)));
   if fid<0
      disp(['Error: Could not open data file ' filenames(ifile,:)]);
      fileerror=1;
      return
   end
   head=fread(fid,4);
   isascii(ifile)= min(head)>=32 & max(head)<=122;

   % Determine filesize by jumping to the end
   fseek(fid,0,'eof');
   filesize(ifile)=ftell(fid);

   % Determine picture size by reading ndim, neqpar, nw, and nx from header
   if isascii(ifile)
      fseek(fid,0,'bof');
      fgetl(fid);
      tmp=sscanf(fgetl(fid),'%d %f %d %d %d',5)';
      ndim=abs(tmp(3)); neqpar=tmp(4); nw=tmp(5);
      nx=sscanf(fgetl(fid),'%d',ndim)'; nxs=prod(nx);

      pictsize(ifile)= 1+79 + 1+7+13+9 + 1+ndim*4 + 1+neqpar*13 + 1+79; %header
      pictsize(ifile)=pictsize(ifile)+ (18*(ndim+nw)+1)*nxs;           %body
   else
      fseek(fid,8+79+4+4+8,'bof');
      ndim=fread(fid,1,'int32'); ndim=abs(ndim);
      neqpar=fread(fid,1,'int32'); nw=fread(fid,1,'int32');
      fread(fid,4); fread(fid,4);
      nx=fread(fid,ndim,'int32'); nxs=prod(nx);

      pictsize(ifile)= 8+79 + 8+4*4+8  + 8+ndim*4 + 8+neqpar*8  + 8+79; %header
      pictsize(ifile)=pictsize(ifile)+ 8*(1+nw+(ndim+nw)*nxs);          %body
   end
   fclose(fid);
end
