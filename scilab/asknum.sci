function [result] = asknum(%prompt,var,nn)

// Ouput variables initialisation (not found in input variables)
result=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Read or display an integer array of nn (optional input param) elements. 

global("Doask");
undef = isempty(var);
if undef then
  result = input(%prompt+"? ");
else
  result = [];
  if Doask then
    // !! L.11: Matlab function sprintf not yet converted, original calling sequence used.
    tmp = sprintf("%g ",var);
    if max(size(mtlb_double(var)))>1 then tmp = "["+trim(tmp)+"] ";end;
    result = input(%prompt+"="+tmp+"? ");
  end;
  if isempty(result) then result = var;end;
end;
// extend or shrunk array if nn is given
if %nargin==3 then
  n = max(size(mtlb_double(result)));
  if mtlb_logic(n,"<",nn) then
    arr = result;
    // ! L.22: real(mtlb_double(nn)) may be replaced by:
    // !    --> mtlb_double(nn) if mtlb_double(nn) is Real.
    result = zeros(1,real(mtlb_double(nn)));
    result = mtlb_i(result,1:n,arr);
    for i = mtlb_imp(n+1,mtlb_double(nn))
      result = mtlb_i(result,i,mtlb_e(arr,n));
    end;
  elseif mtlb_logic(n,">",nn) then
    result = mtlb_e(result,mtlb_imp(1,mtlb_double(nn)));
  end;
end;
if ~undef & ~mtlb_double(Doask) then
  // !! L.32: Matlab function sprintf not yet converted, original calling sequence used.
  tmp = sprintf("%g ",result);
  if max(size(mtlb_double(result)))>1 then tmp = "["+trim(tmp)+"] ";end;
  disp(%prompt+" = "+tmp);
end;
endfunction
