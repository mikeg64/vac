function [arr,n] = str2arr(str,nn)
//function [arr,n] = str2arr(str)

// Ouput variables initialisation (not found in input variables)
arr=[];
n=[];

// Number of arguments in function call
[%nargout,%nargin] = argn(0)

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

// Form a string array from str. In str elements are separated by spaces.
//    arr=str2arr(str)
//    [arr,n]=str2arr(str)
//    arr=str2arr(str,nn)
// The optional nn argument extends the string array arr to have nn elements, 
// and fills up the rest of arr with the last element. The optional n output 
// argument returns the number of elements defined by str. 

l = max(size(mtlb_double(str)));
n = 0;

i = 1;
while i<=l
  while mtlb_isspace(mtlb_e(str,i))
    i = i+1;
    if i>l then break,end;
  end;
  if i>l then break,end;
  i0 = i;
  while ~mtlb_isspace(mtlb_e(str,i))
    i = i+1;
    if i>l then break,end;
  end;
  if n==0 then
    arr = mtlb_e(str,i0:i-1);
  else
    arr = str2mat(arr,mtlb_e(str,i0:i-1));
  end;
  n = n+1;
end;
if n==0 then
  arr = " ";
  n = 1;
end;
if %nargin==2 then
  //nn=argn(2);
  for i = mtlb_imp(n+1,mtlb_double(nn))
    arr = str2mat(arr,arr(n,:));
  end;
  if mtlb_logic(nn,"<",n) then
    // !! L.41: string output can be different from Matlab num2str output.
    disp("Warning: more than "+string(nn)+" values defined by string: "+str);
  end;
end;
endfunction
