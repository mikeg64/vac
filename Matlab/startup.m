clc

disp('Welcome to the MATLAB visualization for VAC');
disp('*******************************************');
disp('By Gabor Toth, October 1996');
disp('********** COMMANDS ***********');
disp('getpict     read snapshots from 1 or more files');
disp('plotfunc    plot functions of last data read by getpict or animate');
disp('animate     read and plot sequence of pictures from 1 or more files');
disp('playmovie   play the movie stored in Movie by animate');
disp('getlog      read the data from a log file');
disp('polargrid   transform coordinates and vector variables to polar');
disp('defaults    clear data, set variables back to defaults values');
disp('********** FUNCTIONS **********');
disp('gradx(f,x)  row-wise    gradient of f');
disp('grady(f,y)  column-wise gradient of f');
disp('********** CUTTING ************');
disp('cut=''2:9,5'' for all functions f plot f(2:9,5) only');
disp('cut=''''      plot the whole function again');
disp('********** ADDING/CHANGING ****');
disp('Doask=1     for confirmation by RETURN, default is Doask=0');
disp('Put extra function definitions into Matlab/get_func.m');
disp('Change the default values in Matlab/defaults.m if you want to');
disp(' ');
defaults;
