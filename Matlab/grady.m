function result=grady(arr,xx)
% Calculate gradient in the second coordinate using central differences
% Copy the gradient at coloumns 2 and n-1 to columns 1 and n
[m,n]=size(arr);
if n<3
   disp('Cannot take X gradient of matrix with less than 3 coloumns');
   return
end
result=zeros(m,n);
result(:,2:n-1)=(arr(:,1:n-2)-arr(:,3:n))./(xx(:,1:n-2)-xx(:,3:n));
result(:,1)=result(:,2);
result(:,n)=result(:,n-1);
