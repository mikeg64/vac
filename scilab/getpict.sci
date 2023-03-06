
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Read the npict-th picture from 1 or more files

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
if mtlb_logic(max(%v02,firstnonsingleton(%v02)),">",1) then
  // ! L.14: mtlb(npict) can be replaced by npict() or npict whether npict is an M-file or not.
  npict = asknum("npict(s) (eg. 2 or [3 4 5])",mtlb(npict),nfile);
else
  npict = npictinfile;
  //dispnum("npict         =",npict);
end;

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
  exec('get_head.sci');
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
  exec('get_body.sci');
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
