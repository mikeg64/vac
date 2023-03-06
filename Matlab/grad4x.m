function result=grad4x(arr,xx)
% Calculate gradient in the first coordinate using 
% 4th order central differences.
% Copy the gradient at rows 3 and m-2 to rows 1,2 and m-1,m
[m,n]=size(arr);
if m<5
   disp('Cannot take 4th order X gradient of matrix with less than 5 rows');
   return
end
result=zeros(m,n);
result(3:m-2,:)=(arr(5:m,:)-8*arr(4:m-1,:)+8*arr(2:m-3,:)-arr(1:m-4,:)) ...
               ./(xx(4:m-1,:)-xx(2:m-3,:))/6;
result(1,:)=result(3,:);
result(2,:)=result(3,:);
result(m-1,:)=result(m-2,:);
result(m,:)=result(m-2,:);
