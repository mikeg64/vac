function [%isascii,filesize,pictsize,fileerror] = studyfile(filenames)

// Ouput variables initialisation (not found in input variables)
%isascii=[];
filesize=[];
pictsize=[];
fileerror=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Determine whether the files are ascii (1) or binary (0) based on 
// the first 4 characters. Determine the filesize and the size of one snapshot.
// Fileerror is set to 1 if any of the files could not be opened. 
// Usage:
//    [isascii,filesize,pictsize,fileerror]=studyfile(filenames)

nfile = size(mtlb_double(filenames),1);
isascii = zeros(1,nfile);
filesize = zeros(1,nfile);
pictsize = zeros(1,nfile);
fileerror = 0;

for ifile = 1:nfile
  fid = mtlb_fopen(trim(filenames(ifile,:)),'r');
  if fid<0 then
    disp("Error: Could not open data file "+filenames(ifile,:));
    fileerror = 1;
    return;
  end;
  // L.21: No simple equivalent, so mtlb_fread() is called.
  head = mtlb_fread(fid,4);
  %isascii = mtlb_i(%isascii,ifile,min(head)>=32 & max(head)<=122);

  // Determine filesize by jumping to the end
  mseek(0,fid,"end");
  filesize = mtlb_i(filesize,ifile,mtell(fid));

  // Determine picture size by reading ndim, neqpar, nw, and nx from header
  if %isascii(ifile) then
    mseek(0,fid,"set");
    %v0_1 = mgetl(fid,1);  if meof()~=0 then %v0_1 = -1;end;  %v0_1;
    %v1_1 = mgetl(fid,1);  if meof()~=0 then %v1_1 = -1;end;  // !! L.32: Matlab function sscanf not yet converted, original calling sequence used.
    tmp = mtlb_t(sscanf(%v1_1,"%d %f %d %d %d",5));
    ndim = abs(mtlb_double(mtlb_e(tmp,3)));  neqpar = mtlb_e(tmp,4);  nw = mtlb_e(tmp,5);
    %v2_1 = mgetl(fid,1);  if meof()~=0 then %v2_1 = -1;end;  // !! L.34: Matlab function sscanf not yet converted, original calling sequence used.
    nx = mtlb_t(sscanf(%v2_1,"%d",ndim));  nxs = mtlb_prod(mtlb_double(nx));
  
    pictsize = mtlb_i(pictsize,ifile,mtlb_a(mtlb_a(mtlb_a(mtlb_a(mtlb_a(1+79+1+7+13+9+1,ndim*4),1),neqpar*13),1),79));  //header
    pictsize = mtlb_i(pictsize,ifile,mtlb_a(pictsize(ifile),mtlb_a(18*mtlb_a(ndim,nw),1)*nxs));  //body
  else
    mseek(8+79+4+4+8,fid,"set");
    // L.40: No simple equivalent, so mtlb_fread() is called.
    ndim = mtlb_fread(fid,1,"int32");  ndim = abs(ndim);
    // L.41: No simple equivalent, so mtlb_fread() is called.
    neqpar = mtlb_fread(fid,1,"int32");  // L.41: No simple equivalent, so mtlb_fread() is called.
    nw = mtlb_fread(fid,1,"int32");
    // L.42: No simple equivalent, so mtlb_fread() is called.
    mtlb_fread(fid,4);  // L.42: No simple equivalent, so mtlb_fread() is called.
    mtlb_fread(fid,4);
    // L.43: No simple equivalent, so mtlb_fread() is called.
    nx = mtlb_fread(fid,ndim,"int32");  nxs = prod(nx,firstnonsingleton(nx));
  
    pictsize = mtlb_i(pictsize,ifile,8+79+8+4*4+8+8+ndim*4+8+neqpar*8+8+79);  //header
    pictsize = mtlb_i(pictsize,ifile,mtlb_a(pictsize(ifile),8*mtlb_a(1+nw,(ndim+nw)*nxs)));  //body
  end;
  mclose(fid);
end;
endfunction
