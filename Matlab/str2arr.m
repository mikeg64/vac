function [arr,n]=str2arr(str,nn)
% Form a string array from str. In str elements are separated by spaces.
%    arr=str2arr(str)
%    [arr,n]=str2arr(str)
%    arr=str2arr(str,nn)
% The optional nn argument extends the string array arr to have nn elements, 
% and fills up the rest of arr with the last element. The optional n output 
% argument returns the number of elements defined by str. 

l=length(str);
n=0;

i=1;
while i<=l
   while isspace(str(i))
     i=i+1;
     if i>l,break,end;
   end
   if i>l,break,end;
   i0=i;
   while ~isspace(str(i))
     i=i+1;
     if i>l,break,end;
   end
   if n==0
      arr=str(i0:i-1);
   else
      arr=str2mat(arr,str(i0:i-1));
   end
   n=n+1;
end
if n==0
   arr=' ';
   n=1;
end
if nargin==2
   for i=n+1:nn
      arr=str2mat(arr,arr(n,:));
   end
   if nn<n
   disp(['Warning: more than ',num2str(nn),' values defined by string: ',str]);
   end
end
