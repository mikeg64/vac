% Plot the variables stored in the memory

disp('======= PLOTTING FUNCTIONS FROM MEMORY ======');
if nfile < 1
   disp('No file has been read yet, run getpict or animate!');
   return
end

if nfile > 1
   disp('More than one files were read...');
   disp(['Probably data is from file ' filenames(nfile,:)]);
   nfile=1;
end
% Restrictions to the last file
ifile=1;
isfirstpict=1;

disp(' ');
disp('======= PLOTTING PARAMETERS =================');

disp(['Variables= ',Variables]);
read_plot_par;

read_limits;
if any(autorange)
   isfirstpict=1;
   get_limits;
   disp(  ['func= ''' func ''' ']);
   dispnum('fmin= ',fmin);
   dispnum('fmax= ',fmax);
end

if ~isempty(multiplot)
   multix=multiplot(1);
   multiy=multiplot(2);
   multidir=multiplot(3);
else
   multix=floor(sqrt(nfunc-1)+1);
   multiy=floor((nfunc-1)/multix+1);
   multidir=0;
end;
npict1=1;

disp('======= PLOT FUNCTIONS ======================');
clf reset; figure(gcf);
ipict=0;
ipict1=0;
iplot=0;
multi0=0;
plot_func;
