function [] = fprintnum(str,arr) 

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Fprint array of integers with a string in front of it, no ending newline

if isempty(arr) then
  tmp = "[]";
else
  // !! L.8: Matlab function sprintf not yet converted, original calling sequence used.
  tmp = sprintf("%g ",arr);
  if max(size(mtlb_double(arr)))>1 then tmp = "["+trim(tmp)+"]";end;
end;
// L.11: No simple equivalent, so mtlb_fprintf() is called.
mtlb_fprintf([str,tmp]);

endfunction
