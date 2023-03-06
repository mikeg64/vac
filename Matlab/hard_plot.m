% Print current figure to printer or to a file while doing animate
% The current value of iplot+1 is added to the filename
% The default filename is Movie/matlab1, Movie/matlab2, ... for iplot=0,1...

if iplot==0
   Printfile=askstr('Printfile (printer/FILENAME) ',Printfile);
   Device=askstr('Device (ps/psc/eps/epsc/gif8...)',Device);
   Orient=askstr('Orient (tall/landscape/portrait)',Orient);
end

if strcmp(Printfile,'printer')
   printfile=''
else
   printfile=[Printfile num2str(iplot+1)];
end

tmp=['print -d' Device ' ' printfile];
fprintf(['Execute: ' tmp ' ?']);
ans=input('Hit RETURN to print, enter 0 to skip, -1 to quit ? ');
if isempty(ans)
   orient(Orient);
   eval(tmp);
elseif ans==-1
   anyerror=1;
end
clear tmp ans printfile;
