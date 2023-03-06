function dispnum(str,arr);
% Fprint array of integers with a string in front of it, no ending newline

if isempty(arr); 
   tmp='[]';
else
   tmp=sprintf('%g ',arr);
   if length(arr)>1; tmp=['[' trim(tmp) ']']; end;
end;
fprintf([str tmp]);

