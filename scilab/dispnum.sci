function [] = dispnum(str,arr) 

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Display array of integers with a string in front of it

if isempty(arr) then
  tmp = "[]";
else
  // !! L.8: Matlab function sprintf not yet converted, original calling sequence used.
  tmp = sprintf("%s ",arr);
  if max(size(mtlb_double(arr)))>1 then tmp = "["+trim(tmp)+"]";end;
end;
disp([str,tmp]);
endfunction
