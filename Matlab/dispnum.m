function dispnum(str,arr);
% Display array of integers with a string in front of it

if isempty(arr); 
   tmp='[]';
else
   tmp=sprintf('%g ',arr);
   if length(arr)>1; tmp=['[' trim(tmp) ']']; end;
end;
disp([str tmp]);
