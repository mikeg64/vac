% Read pictures from one or more VAC data files and plot/animate functions

disp('======= FILE PARAMETERS =====================');
filename=askstr('filename(s)   ',filename);
[filenames,nfile]=str2arr(filename);
[asciifile,filesize,pictsize,fileerror]=studyfile(filenames);
if fileerror; return; end;
dispnum('asciifile(s)  =',asciifile);
npictinfile=floor(filesize./pictsize);
dispnum('npictinfile(s)=',npictinfile);

disp('======= FILE HEADERS ========================');

anygencoord=0;
for ifile=1:nfile
   fid=fopen(trim(filenames(ifile,:)));
   get_head;
   disp(' ');
   disp(['headline=' headline]);
   disp(['Variables='  Variables ' (ndim=' num2str(ndim) ...
                                   ', nw=' num2str(nw) ')']);
   anygencoord=anygencoord | gencoord;
   fclose(fid);
end

disp('======= ANIMATION PARAMETERS ================');
firstpict =asknum('firstpict(s)',firstpict,nfile);
dpict     =asknum('dpict(s)    ',dpict,nfile);
npictmax  =asknum('npictmax    ',npictmax,1);
doanimate =asknum('doanimate   ',doanimate);
dohardplot=asknum('dohardplot  ',dohardplot);
if dohardplot & doanimate
   disp('Since dohardplot=1, doanimate=0 is set to save time');
   doanimate=0;
end

disp('======= PLOTTING PARAMETERS =================');

if anygencoord; read_transform_par; end;
read_plot_par;
read_limits;
if any(autorange)
   for ifile=1:nfile
      fids(ifile)=fopen(trim(filenames(ifile,:)));
   end
   npict=0;
   anyerror=0;
   while npict < npictmax & ~fileerror
      for ifile=1:nfile
         if npict==0
            nextpict=firstpict(ifile);
         else
            nextpict=dpict(ifile);
         end
         fid=fids(ifile);
         fseek(fid,pictsize(ifile)*(nextpict-1),'cof');
         get_head;
         anyerror= fileerror | anyerror;
         if ~anyerror
            get_body;
            isfirstpict= npict==0 & ifile==1;
            get_limits;
         end
      end
      if ~anyerror
         dispnum('ipict: ',firstpict+npict.*dpict);
         npict=npict+1;
      end
   end
   disp(  ['func= ',func]);
   dispnum('fmin= ',fmin);
   dispnum('fmax= ',fmax);
   fclose('all');
else
   disp(' ');
   npict=max(floor((npictinfile-firstpict)./dpict+1));
   if npict>npictmax; npict=npictmax; end;
   disp(['npict = ',num2str(npict)]);
end

disp('======= PLOT FUNCTIONS ======================');

if ~isempty(multiplot)
   multix=multiplot(1);
   multiy=multiplot(2);
   multidir=multiplot(3);
   npict1=floor((multix*multiy)/(nfunc*nfile));
   if npict1==0; npict1=1; end;
elseif nfile==1
   multix=floor(sqrt(nfunc-1)+1);
   multiy=floor((nfunc-1)/multix+1);
   multidir=0;
   npict1=1;
else
   multix=nfile;
   multiy=nfunc;
   multidir=1;
   npict1=1;
end

nplot=floor((npict-1)/npict1)+1;

clf reset; figure(gcf); 

for ifile=1:nfile
   fids(ifile)=fopen(trim(filenames(ifile,:)));
end
ipict=0;
ipict1=0;
iplot=0;
multi0=0;
anyerror=0;
while ipict < npict & ~anyerror
   if ipict1==0; clf; end;
   for ifile=1:nfile
      if ipict==0
         nextpict=firstpict(ifile);
      else
         nextpict=dpict(ifile);
      end
      fid=fids(ifile);
      fseek(fid,pictsize(ifile)*(nextpict-1),'cof');
      get_head;
      anyerror= fileerror | anyerror;
      if ~anyerror
         get_body;
         isfirstpict= ipict==0 & ifile==1;
         plot_func;
      end;
   end;
   if ~anyerror;
      dispnum('ipict: ',firstpict+ipict.*dpict);

      if ipict1==npict1-1 | ipict==npict-1
         if doanimate & npict>npict1
            if iplot==0;Movie=moviein(nplot,gcf);end
            Movie(:,iplot+1)=getframe(gcf);
         elseif ~dohardplot
            fprintf('Press any key to proceed...');pause;fprintf('\n');
         end
         if dohardplot
            hard_plot;
         end
      end
      ipict1=ipict1+1;
      if ipict1 >= npict1; ipict1=0; multi0=0; iplot=iplot+1; end;
      ipict=ipict+1;
   end
end
fclose('all');
if doanimate & npict>npict1
  if ipict<npict;Movie=Movie(:,1:ipict);end;
  playmovie;
end
npict=ipict;
