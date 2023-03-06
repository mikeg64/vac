% Read limits for functions from input for which autorange is false

disp('======= Determine plotting range(s) =========');

if ~all(autorange)
   disp('Minimum and maximum value(s) for function(s) (eg. 0.3 or [0.1 2])');
   disp('(arbitrary values for those with autorange==1 set)');
   disp(      ['func(s) = ''',func ''' ']);
   fmin=asknum('fmin(s)',fmin,nfunc);
   fmax=asknum('fmax(s)',fmax,nfunc);
else
   fmin=zeros(1,nfunc);
   fmax=zeros(1,nfunc);
end
