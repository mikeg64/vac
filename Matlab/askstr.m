function [result]= askstr(prompt,var)

global Doask;

if isempty(var)
   result=input([prompt '? ']);
elseif Doask
   result=input([prompt '=''' var ''' ? ']);
   if isempty(result); result=var; end;
else
   disp([prompt ' = ''' var ''' ']);
   result=var;
end

