function result= asknum(prompt,var,nn)
% Read or display an integer array of nn (optional input param) elements. 

global Doask;
undef= isempty(var);
if undef 
   result=input([prompt,'? ']);
else
   result=[];
   if Doask; 
      tmp=sprintf('%g ',var); 
      if(length(var)>1);tmp=['[',trim(tmp),'] ',];end;
      result=input([prompt,'=',tmp,'? ']); 
   end;
   if isempty(result); result=var; end;
end
% extend or shrunk array if nn is given
if nargin==3
   n=length(result);
   if n<nn
      arr=result;
      result=zeros(1,nn);
      result(1:n)=arr;
      for i=n+1:nn
         result(i)=arr(n);
      end
   elseif n>nn
      result=result(1:nn);
   end
end
if ~undef & ~Doask
   tmp=sprintf('%g ',result); 
   if(length(result)>1);tmp=['[',trim(tmp),'] ',];end;
   disp([prompt ' = ' tmp]);
end
