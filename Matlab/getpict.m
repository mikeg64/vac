% Read the npict-th picture from 1 or more files

filename=askstr('filename(s)  ',filename);
[filenames,nfile]=str2arr(filename);
[asciifile,filesize,pictsize,fileerror]=studyfile(filenames);
if fileerror; return; end;
dispnum('asciifile(s)  =',asciifile);
disp(' ');
npictinfile=floor(filesize./pictsize);
fprintnum('npictinfile(s)=',npictinfile);
if ~isempty(npict); fprintf('\n'); end;

if max(npictinfile)>1
   npict=asknum('npict(s) (eg. 2 or [3 4 5])',npict,nfile);
else
   npict=npictinfile;
   dispnum('npict         =',npict);
end

for ifile=1:nfile
   if npict(ifile)>npictinfile(ifile)
      disp(' ');
      disp(['Reducing npict=' num2str(npict(ifile)) ' for file ' ...
         num2str(ifile) ' to npictinfile=' num2str(npictinfile(ifile)) ' !']);
      npict(ifile)=npictinfile(ifile);
   end;
   fid=fopen(trim(filenames(ifile,:)));
   fseek(fid,pictsize(ifile)*(npict(ifile)-1),'bof');
   get_head;
   disp(' ');
   disp(['headline=''',headline,''' ']);
   if ~isempty(physics); disp(['physics =''',physics,''' ']); end;
   if fileerror
      fprintf('Error: Could not read the %d',dpict);
      disp(['. picture from data file ' filenames(ifile,:)]);
      return
   end
   fprintf('it=%d time=%f ndim=%d neqpar=%d nw=%d\n',it,time,ndim,neqpar,nw);
   dispnum('eqpar=',eqpar);
   dispnum('nx   =',nx);
   if gencoord; read_transform_par; end;
   if nfile>1
      disp(['Reading:',Variables,' (with _',num2str(ifile),')']);
   else
      disp(['Reading:',Variables]);
   end;
   get_body;
   fclose(fid);
   if nfile>1
      % Rename the variables to rho_1,rho_2 etc for more than one file
      for i=1:size(variables,1)
        eval([trim(variables(i,:)) '_' num2str(ifile) '=' variables(i,:) ';']);
      end
   end
end
