function [arr,n]=comma2arr(str,nn)
% Form a string array from str. In str elements are separated by semicolons.
%    arr=comma2arr(str)
%    [arr,n]=comma2arr(str)
%    arr=comma2arr(str,nn)
% The optional nn argument extends the string array arr to length nn, 
% and fills up the rest of arr with the last element. The optional n output 
% contains the number of elements defined by str. Leading blanks are cut off
% from the array elements, thus comma2arr('a; b; c; etc') will form the matrix
% arr=['a  ';'b  ';'c  ';'etc']

l=length(str);
n=0;
sep=';';

i=1;
while i<=l
   while isspace(str(i)) | str(i)==sep
     i=i+1;
     if i>l,break,end;
   end
   if i>l,break,end;
   i0=i;
   while str(i)~=sep
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
