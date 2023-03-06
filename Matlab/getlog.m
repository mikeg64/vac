% Read the ASCII log file

logfilename=askstr('logfilename ',logfilename);
fid=fopen(logfilename);
if fid<0
   disp(['Error: Could not open log file ',logfilename]);
   return
end

headline=trim(fgetl(fid));
disp(['headline     = ',headline]);
wlogname=trim(fgetl(fid));
[wlognames,nwlog]=str2arr(wlogname);

disp('Reading the wlog matrix...');
wlog=fscanf(fid,'%f',[nwlog,Inf])';

fprintf('Number of rows in wlog: %d\n',size(wlog,1));
disp('Extracting columns of wlog into arrays:');
disp(wlogname);
for iwlog=1:nwlog
  eval([wlognames(iwlog,:),'=wlog(:,iwlog);']);
end
