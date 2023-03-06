function []=mygetpictfunc()
// Display mode
//mode(0);

// Display warning for floating point exception
//ieee(1);
  labindex=1;
  numlabs=1;
  nr=100;
  ne=400+8;
  nnepp=ne/numlabs;
  startne=400+(labindex-1)*nnepp;
  
  if labindex==numlabs
    finishne=ne;  
  else
    finishne=startne+nnepp;
  end
// Read the npict-th picture from 1 or more files
% File parameters
filename='';
physics='';
phys='';
ndir=[];
logfilename='';

% Transformation parameters
Transform='';
nxreg=[];
nxregold=[];
polar_r='';
polar_phi='';

// Function parameters
func='';
autorange=[];
fmin=[];
fmax=[];

// Plotting parameters
cut='';
plotmode='';
plottitle='default';
View=[-37.5 30];
Colorbar=0;
Shading='flat';
Contourlevel=30;
Contourstyle='g-';
Quiverscale=1;

// multiplot=[] gives the default number of subplots depending on nfile,nfunc
// multiplot=[3,2] defines 3 by 2 subplots
multiplot=[];

// Number of info items on bottom and in the header
Bottomline=2;
Headerline=2;

// Animation parameters
npict=[];
firstpict=1;
dpict=1;
npictmax=100;
doanimate=1;

// Printing parameters
dohardplot=0;
Printfile='Movie/matlab';
Device='ps';
Orient='landscape';

//mypars

// ! L.3: mtlb(filename) can be replaced by filename() or filename whether filename is an M-file or not.
filename='../data/grav2.out';
//filename = askstr("filename(s)  ",mtlb(filename));
[filenames,nfile] = str2arr(filename);
[asciifile,filesize,pictsize,fileerror] = studyfile(filenames);
if fileerror then return;end;
//dispnum("asciifile(s)  =",asciifile);
disp(" ");
npictinfile = floor(filesize ./pictsize);
// !! L.10: Unknown function fprintnum not converted, original calling sequence used.
//fprintnum("npictinfile(s)=",npictinfile);
// ! L.11: mtlb(npict) can be replaced by npict() or npict whether npict is an M-file or not.

if ~isempty(mtlb(npict)) then // L.11: No simple equivalent, so mtlb_fprintf() is called.
 mtlb_fprintf("\n");end;

%v02 = npictinfile;
//if mtlb_logic(max(%v02,firstnonsingleton(%v02)),">",1) then
  // ! L.14: mtlb(npict) can be replaced by npict() or npict whether npict is an M-file or not.
//  npict = asknum("npict(s) (eg. 2 or [3 4 5])",mtlb(npict),nfile);
//else
//  npict = npictinfile;
  //dispnum("npict         =",npict);
//end;

for npict=startne:finishne

%npict=415;
outname=sprintf('../data/ascdata/testascmat_%d.out',npict);
imoutname=sprintf('../data/plots/testascmat_%d.jpg',npict);   



for ifile = 1:nfile
  if mtlb_logic(mtlb_e(npict,ifile),">",npictinfile(ifile)) then
    disp(" ");
  
    // !! L.24: string output can be different from Matlab num2str output.
    // !! L.24: string output can be different from Matlab num2str output.
    // !! L.24: string output can be different from Matlab num2str output.
    disp("Reducing npict="+string(mtlb_e(npict,ifile))+" for file "+string(ifile)+" to npictinfile="+string(npictinfile(ifile))+" !");
    npict = mtlb_i(npict,ifile,npictinfile(ifile));
  end;
  fid = mtlb_fopen(trim(filenames(ifile,:)),'r');
  mseek(pictsize(ifile)*mtlb_s(mtlb_e(npict,ifile),1),fid,"set");
  // !! L.29: Unknown function get_head not converted, original calling sequence used.
  //exec('get_head.sci');
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
  
  
  
  
  
  
  
  disp(" ");
  // ! L.31: mtlb(headline) can be replaced by headline() or headline whether headline is an M-file or not.
  disp("headline=''"+mtlb(headline)+"'' ");
  // ! L.32: mtlb(physics) can be replaced by physics() or physics whether physics is an M-file or not.

  if ~isempty(mtlb(physics)) then // ! L.32: mtlb(physics) can be replaced by physics() or physics whether physics is an M-file or not.
   disp("physics =''"+mtlb(physics)+"'' ");end;
  if fileerror then
    // ! L.34: mtlb(dpict) can be replaced by dpict() or dpict whether dpict is an M-file or not.
    // L.34: No simple equivalent, so mtlb_fprintf() is called.
    mtlb_fprintf("Error: Could not read the %d",mtlb(dpict));
    disp(". picture from data file "+filenames(ifile,:));
    return;
  end;
  // ! L.38: mtlb(it) can be replaced by it() or it whether it is an M-file or not.
  // ! L.38: mtlb(time) can be replaced by time() or time whether time is an M-file or not.
  // ! L.38: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
  // ! L.38: mtlb(neqpar) can be replaced by neqpar() or neqpar whether neqpar is an M-file or not.
  // ! L.38: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.
  // L.38: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf("it=%d time=%f ndim=%d neqpar=%d nw=%d\n",mtlb(it),mtlb(time),mtlb(ndim),mtlb(neqpar),mtlb(nw));
  // ! L.39: mtlb(eqpar) can be replaced by eqpar() or eqpar whether eqpar is an M-file or not.
  //dispnum("eqpar=",mtlb(eqpar));
  // ! L.40: mtlb(nx) can be replaced by nx() or nx whether nx is an M-file or not.
  //dispnum("nx   =",mtlb(nx));
  // ! L.41: mtlb(gencoord) can be replaced by gencoord() or gencoord whether gencoord is an M-file or not.

  if mtlb(gencoord) then // !! L.41: Unknown function read_transform_par not converted, original calling sequence used.
   read_transform_par;end;
  if nfile>1 then
    // ! L.43: mtlb(Variables) can be replaced by Variables() or Variables whether Variables is an M-file or not.
    // !! L.43: string output can be different from Matlab num2str output.
    disp("Reading:"+mtlb(Variables)+" (with _"+string(ifile)+")");
  else
    // ! L.45: mtlb(Variables) can be replaced by Variables() or Variables whether Variables is an M-file or not.
    disp("Reading:"+mtlb(Variables));
  end;
  // !! L.47: Unknown function get_body not converted, original calling sequence used.
  //exec('get_body.sci');
  
  // Read a snapshot from a VAC data file (binary or ascii).
// First X and w are read, (X is capitalized to avoid possible conflict with 
// the name of its first component). Transform X and w according to Transform
// if generalized coordinates are found.
// Extract the variables given in the xwname string array read by get_head.

// ! L.7: mtlb(nxs) can be replaced by nxs() or nxs whether nxs is an M-file or not.
// ! L.7: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
// ! L.7: real(mtlb_double(mtlb(nxs))) may be replaced by:
// !    --> mtlb_double(mtlb(nxs)) if mtlb_double(mtlb(nxs)) is Real.
// ! L.7: real(mtlb_double(mtlb(ndim))) may be replaced by:
// !    --> mtlb_double(mtlb(ndim)) if mtlb_double(mtlb(ndim)) is Real.
X = zeros(real(mtlb_double(mtlb(nxs))),real(mtlb_double(mtlb(ndim))));
// ! L.8: mtlb(nxs) can be replaced by nxs() or nxs whether nxs is an M-file or not.
// ! L.8: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.
// ! L.8: real(mtlb_double(mtlb(nxs))) may be replaced by:
// !    --> mtlb_double(mtlb(nxs)) if mtlb_double(mtlb(nxs)) is Real.
// ! L.8: real(mtlb_double(mtlb(nw))) may be replaced by:
// !    --> mtlb_double(mtlb(nw)) if mtlb_double(mtlb(nw)) is Real.
w = zeros(real(mtlb_double(mtlb(nxs))),real(mtlb_double(mtlb(nw))));
// ! L.9: mtlb(asciifile) can be replaced by asciifile() or asciifile whether asciifile is an M-file or not.
// ! L.9: mtlb(ifile) can be replaced by ifile() or ifile whether ifile is an M-file or not.
// !! L.9: Unknown function asciifile not converted, original calling sequence used.

if asciifile(mtlb(ifile)) then
  // ! L.10: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // ! L.10: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
  // ! L.10: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.
  // ! L.10: mtlb(nxs) can be replaced by nxs() or nxs whether nxs is an M-file or not.
  // L.10: No simple equivalent, so mtlb_fscanf() is called.
  xw = mtlb_fscanf(mtlb(fid),"%f",[mtlb_a(mtlb(ndim),mtlb(nw)),mtlb(nxs)]);
  // ! L.11: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
  X = (xw(mtlb_imp(1,mtlb_double(mtlb(ndim))),:))';
  // ! L.12: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
  // ! L.12: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.
  // ! L.12: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.
  w = (xw(mtlb_imp(mtlb_a(mtlb(ndim),1),mtlb_a(mtlb(nw),mtlb(ndim))),:))';
  clear("xw");
  // ! L.14: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  %v0_1 = mgetl(mtlb(fid),1);  if meof()~=0 then %v0_1 = -1;end;  %v0_1;
else
  // ! L.16: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.16: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.17: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

  for idim = mtlb_imp(1,mtlb_double(mtlb(ndim)))
    // ! L.18: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
    // ! L.18: mtlb(nxs) can be replaced by nxs() or nxs whether nxs is an M-file or not.
    // L.18: No simple equivalent, so mtlb_fread() is called.
    X(:,idim) = mtlb_fread(mtlb(fid),mtlb(nxs),"float64");
  end;
  // ! L.20: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
  // L.20: No simple equivalent, so mtlb_fread() is called.
  mtlb_fread(mtlb(fid),4);
  // ! L.21: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.

  for iw = mtlb_imp(1,mtlb_double(mtlb(nw)))
    // ! L.22: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
    // L.22: No simple equivalent, so mtlb_fread() is called.
    mtlb_fread(mtlb(fid),4);
    // ! L.23: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
    // ! L.23: mtlb(nxs) can be replaced by nxs() or nxs whether nxs is an M-file or not.
    // L.23: No simple equivalent, so mtlb_fread() is called.
    w(:,iw) = mtlb_fread(mtlb(fid),mtlb(nxs),"float64");
    // ! L.24: mtlb(fid) can be replaced by fid() or fid whether fid is an M-file or not.
    // L.24: No simple equivalent, so mtlb_fread() is called.
    mtlb_fread(mtlb(fid),4);
  end;
end;

// extract variables from X into variables named after the strings in xnames
// ! L.29: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

for idim = mtlb_imp(1,mtlb_double(mtlb(ndim)))
  // ! L.30: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

  if mtlb_logic(mtlb(ndim),"==",2) then
    // ! L.31: mtlb(nx1) can be replaced by nx1() or nx1 whether nx1 is an M-file or not.
    // ! L.31: mtlb(nx2) can be replaced by nx2() or nx2 whether nx2 is an M-file or not.
    // !! L.31: WARNING: Matlab reshape() suppresses singleton higher dimension, it is not the case for matrix.
    tmp = matrix(X(:,idim),mtlb(nx1),mtlb(nx2));
  else
    tmp = X(:,idim);
  end;
  // !! L.35: Unknown function xnames not converted, original calling sequence used.
  disp("this bit  1");
  //mtlb_eval(xnames(idim,:)+"=tmp;");
  sprintf("%s%s",xnames(idim,:),"=tmp;");
end;
disp("this bit  2");
// ! L.39: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

if mtlb_logic(mtlb(ndim),"==",1) then
  xx = X;
elseif mtlb_logic(mtlb(ndim),"==",2) then
  // ! L.40: mtlb(nx1) can be replaced by nx1() or nx1 whether nx1 is an M-file or not.
  // ! L.40: mtlb(nx2) can be replaced by nx2() or nx2 whether nx2 is an M-file or not.
  // !! L.40: WARNING: Matlab reshape() suppresses singleton higher dimension, it is not the case for matrix.
  xx = matrix(X(:,1),mtlb(nx1),mtlb(nx2));
  // ! L.41: mtlb(nx1) can be replaced by nx1() or nx1 whether nx1 is an M-file or not.
  // ! L.41: mtlb(nx2) can be replaced by nx2() or nx2 whether nx2 is an M-file or not.
  // !! L.41: WARNING: Matlab reshape() suppresses singleton higher dimension, it is not the case for matrix.
  yy = matrix(X(:,2),mtlb(nx1),mtlb(nx2));
end;

// extract variables from w into variables named after the strings in wnames
// ! L.45: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.

for iw = mtlb_imp(1,mtlb_double(mtlb(nw)))
  // ! L.46: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

  if mtlb_logic(mtlb(ndim),"==",2) then
    // ! L.47: mtlb(nx1) can be replaced by nx1() or nx1 whether nx1 is an M-file or not.
    // ! L.47: mtlb(nx2) can be replaced by nx2() or nx2 whether nx2 is an M-file or not.
    // !! L.47: WARNING: Matlab reshape() suppresses singleton higher dimension, it is not the case for matrix.
    tmp = matrix(w(:,iw),mtlb(nx1),mtlb(nx2));
  else
    tmp = w(:,iw);
  end;
  // !! L.51: Unknown function wnames not converted, original calling sequence used.
  //mtlb_eval(wnames(iw,:)+"=tmp;");
  sprintf("%s%s",wnames(iw,:),"=tmp;")
end;
clear("tmp");

// ! L.55: mtlb(gencoord) can be replaced by gencoord() or gencoord whether gencoord is an M-file or not.
// ! L.55: mtlb(ndim) can be replaced by ndim() or ndim whether ndim is an M-file or not.

if mtlb_double(mtlb(gencoord)) & mtlb_logic(mtlb(ndim),"==",2) then
  // ! L.58: mtlb(Transform) can be replaced by Transform() or Transform whether Transform is an M-file or not.

  if mtlb_strcmp(mtlb(Transform),"polar") then
    // !! L.57: Unknown function polargrid not converted, original calling sequence used.
    polargrid;
  elseif mtlb_strcmp(mtlb(Transform),"regular") then
    // !! L.59: Unknown function regulargrid not converted, original calling sequence used.
    regulargrid;
  end;
end;

  
  
  
  
  mclose(fid);
  if nfile>1 then
    // Rename the variables to rho_1,rho_2 etc for more than one file
    // ! L.51: mtlb(variables) can be replaced by variables() or variables whether variables is an M-file or not.
  
    for i = 1:size(mtlb_double(mtlb(variables)),1)
      // ! L.52: mtlb(variables) can be replaced by variables() or variables whether variables is an M-file or not.
      // !! L.52: Unknown function variables not converted, original calling sequence used.
      // !! L.52: string output can be different from Matlab num2str output.
      // ! L.52: mtlb(variables) can be replaced by variables() or variables whether variables is an M-file or not.
      // !! L.52: Unknown function variables not converted, original calling sequence used.
      mtlb_eval(trim(variables(i,":"))+"_"+string(ifile)+"="+variables(i,":")+";");
    end;
  end;
end;


//reduce the data and save in ascii format
[nr,nc]=size(x);
ntot=nr*nc;

ascdat=zeros(nr,nc,12);
ascdat(:,:,1)=x(:,:);
ascdat(:,:,2)=y(:,:);
ascdat(:,:,3)=h(:,:);
ascdat(:,:,4)=m1(:,:);
ascdat(:,:,5)=m2(:,:);
ascdat(:,:,6)=e(:,:);
ascdat(:,:,7)=b1(:,:);
ascdat(:,:,8)=b2(:,:);
ascdat(:,:,9)=eb(:,:);
ascdat(:,:,10)=rhob(:,:);
ascdat(:,:,11)=bg1(:,:);
ascdat(:,:,12)=bg2(:,:);




end; //looping over pictures

endfunction

