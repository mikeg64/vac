
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Read the transformation parameters for irregular grids

// ! L.3: mtlb(Transform) can be replaced by Transform() or Transform whether Transform is an M-file or not.
Transform = askstr("Transform (regular/polar/none)",mtlb(Transform));
if mtlb_strcmp(Transform,"regular") then
  // ! L.5: mtlb(nxreg) can be replaced by nxreg() or nxreg whether nxreg is an M-file or not.
  nxreg = asknum("nxreg (eg. 100 or [50 60])",mtlb(nxreg),2);
elseif mtlb_strcmp(Transform,"polar") then
  disp("Radial and phi component(s) of vector variables among");
  // L.8: No simple equivalent, so mtlb_fprintf() is called.
  mtlb_fprintf("wnames                           = ");
  // ! L.9: mtlb(nw) can be replaced by nw() or nw whether nw is an M-file or not.

  for iw = mtlb_imp(1,mtlb_double(mtlb(nw))) // !! L.9: Unknown function wnames not converted, original calling sequence used.
   // L.9: No simple equivalent, so mtlb_fprintf() is called.
   mtlb_fprintf(trim(wnames(iw,":"))+" ");end;
  // ! L.10: mtlb(polar_r) can be replaced by polar_r() or polar_r whether polar_r is an M-file or not.

  if ~isempty(mtlb(polar_r)) then // L.10: No simple equivalent, so mtlb_fprintf() is called.
   mtlb_fprintf("\n");end;
  // ! L.11: mtlb(polar_r) can be replaced by polar_r() or polar_r whether polar_r is an M-file or not.
  polar_r = askstr("polar_r(s)   (eg. mr br     ...)",mtlb(polar_r));
  // ! L.12: mtlb(polar_phi) can be replaced by polar_phi() or polar_phi whether polar_phi is an M-file or not.
  polar_phi = askstr("polar_phi(s) (eg. mphi bphi ...)",mtlb(polar_phi));
  polar_rs = str2arr(polar_r);
  polar_phis = str2arr(polar_phi);
elseif ~mtlb_strcmp(Transform,"none") then
  disp("Incorrect value "+Transform+" is replaced by ''none''");
  Transform = "none"
end;
