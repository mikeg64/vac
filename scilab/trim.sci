function [result] = trim(str)

// Ouput variables initialisation (not found in input variables)
result=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//Trim trailing blanks from str

i = max(size(mtlb_double(str)));
while mtlb_isspace(mtlb_e(str,i))
  i = i-1;
end;
result = mtlb_e(str,1:i);
endfunction
