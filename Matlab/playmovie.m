% Plays movie saved into Movie by animate.

if isempty(Movie)
   disp('Movie is empty, run animate first!');
   return;
end;

disp('======= PLAYING THE MOVIE ===================');
disp('How many times to play the movie? Optionally give speed as 2 decimals,');
disp('e.g. 5.04 for 5 times with 4 frames/sec. The default speed is 12.');

while 1
   tmp=[];
   while isempty(tmp)
      tmp=input('nmovie (0 to quit, -1 for frame by frame)? ');
   end
   nmovie=floor(tmp);
   if nmovie==0;break;end;
   clf reset;figure(gcf);
   if nmovie<0
      disp(['   nframe=',num2str(size(Movie,2))]);
      iframe=1;
      while iframe>0
         movie(gcf,Movie(:,iframe));
         tmp=input(['   iframe=',num2str(iframe), ...
                    ' (new value/return for next/0 to quit)? ']);
         if isempty(tmp);iframe=iframe+1;else;iframe=tmp;end;
         if iframe>size(Movie,2);iframe=size(Movie,2);end;
      end
      clear iframe;
   else
      speed=100*(tmp-nmovie);
      if speed<0.1;speed=12.5;end;
      movie(gcf,Movie,nmovie,speed);
   end
end

clear tmp nmovie speed;
