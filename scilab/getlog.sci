
// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Read the ASCII log file

// ! L.3: mtlb(logfilename) can be replaced by logfilename() or logfilename whether logfilename is an M-file or not.
logfilename = askstr("logfilename ",mtlb(logfilename));
fid = mtlb_fopen(logfilename);
if fid<0 then
  disp("Error: Could not open log file "+logfilename);
  return;
end;

%v0 = mgetl(fid,1);if meof()~=0 then %v0 = -1;end;headline = trim(%v0);
disp("headline     = "+headline);
%v0 = mgetl(fid,1);if meof()~=0 then %v0 = -1;end;wlogname = trim(%v0);
[wlognames,nwlog] = str2arr(wlogname);

disp("Reading the wlog matrix...");
// L.16: No simple equivalent, so mtlb_fscanf() is called.
wlog = mtlb_fscanf(fid,"%f",[nwlog,%inf])';

// L.18: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf("Number of rows in wlog: %d\n",size(wlog,1));
disp("Extracting columns of wlog into arrays:");
disp(wlogname);
for iwlog = 1:nwlog
  mtlb_eval(wlognames(iwlog,:)+"=wlog(:,iwlog);");
end;
