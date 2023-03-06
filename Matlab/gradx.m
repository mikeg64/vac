function result=gradx(arr,xx)
% Calculate gradient in the first coordinate using central differences
% Take the gradient at the 2-nd and n
[m,n]=size(arr);
if m<3
   disp('Cannot take X gradient of matrix with less than 3 rows');
   return
end
result=zeros(m,n);
result(2:m-1,:)=(arr(1:m-2,:)-arr(3:m,:))./(xx(1:m-2,:)-xx(3:m,:));
result(1,:)=result(2,:);
result(m,:)=result(m-1,:);
