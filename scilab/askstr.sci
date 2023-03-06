function [result] = askstr(%prompt,var)

// Ouput variables initialisation (not found in input variables)
result=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);


global("Doask");

if isempty(var) then
  result = input(%prompt+"? ");
elseif Doask then
  result = input(%prompt+"=''"+var+"'' ? ");
  if isempty(result) then result = var;end;
else
  disp(%prompt+" = ''"+var+"'' ");
  result = var;
end;

endfunction
