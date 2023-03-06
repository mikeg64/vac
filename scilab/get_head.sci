
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Read header of a VAC data file

fileerror = 0;
// ! L.4: mtlb(asciifile) can be replaced by asciifile() or asciifile whether asciifile is an M-file or not.
// ! L.4: mtlb(ifile) can be replaced by ifile() or ifile whether ifile is an M-file or not.
// !! L.4: Unknown function asciifile not converted, original calling sequence used.

if asciifile(mtlb(ifile)) then
  // ! L.5: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v0_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v0_1 = -1;end;  headline = trim(%v0_1);
  if ~type(headline)==10 then fileerror = 1; return;end;
  // ! L.7: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v1_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v1_1 = -1;end;  // !! L.7: Matlab function sscanf not yet converted, original calling sequence used.
  tmp = mtlb_t(sscanf(%v1_1,"%d %f %d %d %d",5));
  it = mtlb_e(tmp,1);  time = mtlb_e(tmp,2);  ndim = mtlb_e(tmp,3);  neqpar = mtlb_e(tmp,4);  nw = mtlb_e(tmp,5);  clear("tmp");
  gencoord = mtlb_logic(ndim,"<",0);  ndim = abs(mtlb_double(ndim));
  // ! L.10: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v2_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v2_1 = -1;end;  // !! L.10: Matlab function sscanf not yet converted, original calling sequence used.
  nx = mtlb_t(sscanf(%v2_1,"%d",ndim));
  // ! L.11: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v3_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v3_1 = -1;end;  // !! L.11: Matlab function sscanf not yet converted, original calling sequence used.
  eqpar = mtlb_t(sscanf(%v3_1,"%f",neqpar));
  // ! L.12: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v4_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v4_1 = -1;end;  Variables = trim(%v4_1);
else
  // ! L.14: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.14: No simple equivalent, so mtlb_fread() is called.
  [tmp,ntmp] = mtlb_fread(mtlb(fid),4);
  if ntmp<4 then fileerror = 1; return;end;
  // ! L.16: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.16: No simple equivalent, so mtlb_fread() is called.
  headline = trim(ascii(mtlb_fread(mtlb(fid),79)'));
  // ! L.17: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.17: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);

  // ! L.19: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.19: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.20: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.20: No simple equivalent, so mtlb_fread() is called.
  it = mtlb_fread(mtlb(fid),1,"int32");  // ! L.20: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.20: No simple equivalent, so mtlb_fread() is called.
  time = mtlb_fread(mtlb(fid),1,"float64");
  // ! L.21: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.21: No simple equivalent, so mtlb_fread() is called.
  ndim = mtlb_fread(mtlb(fid),1,"int32");
  gencoord = ndim<0;  ndim = abs(ndim);
  // ! L.23: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.23: No simple equivalent, so mtlb_fread() is called.
  neqpar = mtlb_fread(mtlb(fid),1,"int32");  // ! L.23: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.23: No simple equivalent, so mtlb_fread() is called.
  nw = mtlb_fread(mtlb(fid),1,"int32");
  // ! L.24: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.24: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);

  // ! L.26: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.26: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.27: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.27: No simple equivalent, so mtlb_fread() is called.
  nx = mtlb_fread(mtlb(fid),ndim,"int32");
  // ! L.28: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.28: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);

  // ! L.30: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.30: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.31: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.31: No simple equivalent, so mtlb_fread() is called.
  eqpar = mtlb_fread(mtlb(fid),neqpar,"float64");
  // ! L.32: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.32: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);

  // ! L.34: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.34: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.35: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.35: No simple equivalent, so mtlb_fread() is called.
  Variables = trim(ascii(mtlb_fread(mtlb(fid),79)'));
  // ! L.36: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.36: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
end;

// Extract physics from headline if it is defined
i = max(size(mtlb_double(headline)));
while i>1 & mtlb_logic(mtlb_e(headline,i),"~=","_") i = i-1;end;
if mtlb_logic(mtlb_e(headline,i),"==","_") then
  physics = mtlb_e(headline,i+1:max(size(mtlb_double(headline))));
  headline = mtlb_e(headline,1:i-1);
end;
// Extraxt number of vector components NDIR from last character of physics
// and extract phys without the number of dimesions and components
if ~isempty(physics) then
  ndir = evstr(mtlb_e(physics,max(size(mtlb_double(physics)))));
  phys = mtlb_e(physics,1:max(size(mtlb_double(physics)))-2);
end;

// Extract info from names of Variables
[variables,ntmp] = str2arr(Variables);
xnames = variables(mtlb_imp(1,ndim),:);
wnames = variables(mtlb_imp(mtlb_a(ndim,1),mtlb_a(ndim,nw)),:);

// It can optionally contain the names of the equation parameters
if mtlb_logic(ntmp,"==",mtlb_a(mtlb_a(ndim,nw),neqpar)) then
  eqparnames = variables(mtlb_imp(mtlb_a(mtlb_a(ndim,nw),1),ntmp),:);
  for ieqpar = mtlb_imp(1,mtlb_double(neqpar))
    mtlb_eval(eqparnames(ieqpar,:)+"= eqpar(ieqpar);");
  end;
end;
clear("ntmp");

nxs = mtlb_prod(mtlb_double(nx));
if mtlb_logic(ndim,"==",1) then
  nx1 = nx;
  nx2 = 1;
  nx3 = 1;
end;
if mtlb_logic(ndim,"==",2) then
  nx1 = mtlb_e(nx,1);
  nx2 = mtlb_e(nx,2);
  nx3 = 1;
end;
if mtlb_logic(ndim,"==",3) then
  nx1 = mtlb_e(nx,1);
  nx2 = mtlb_e(nx,2);
  nx3 = mtlb_e(nx,3);
end;
