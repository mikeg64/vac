
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

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
