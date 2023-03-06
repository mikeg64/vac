% Default settings for all the macros

clear;

% Confirmation of parameters already set
global Doask;
Doask=0;

% File parameters
filename='';
physics='';
phys='';
ndir=[];
logfilename='';

% Transformation parameters
Transform='';
nxreg=[];
nxregold=[];
polar_r='';
polar_phi='';

% Function parameters
func='';
autorange=[];
fmin=[];
fmax=[];

% Plotting parameters
cut='';
plotmode='';
plottitle='default';
View=[-37.5 30];
Colorbar=0;
Shading='flat';
Contourlevel=30;
Contourstyle='g-';
Quiverscale=1;

% multiplot=[] gives the default number of subplots depending on nfile,nfunc
% multiplot=[3,2] defines 3 by 2 subplots
multiplot=[];

% Number of info items on bottom and in the header
Bottomline=2;
Headerline=2;

% Animation parameters
npict=[];
firstpict=1;
dpict=1;
npictmax=100;
doanimate=1;

% Printing parameters
dohardplot=0;
Printfile='Movie/matlab';
Device='ps';
Orient='landscape';
