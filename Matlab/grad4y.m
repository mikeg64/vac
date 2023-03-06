function result=grad4y(arr,yy)
% Calculate gradient in the second coordinate using 
% 4th order central differences. 
% Copy the gradient at coloumns 3 and n-2 to columns 1,2 and n-1,n
[m,n]=size(arr);
if n<5
 disp('Cannot take 4th order Y gradient of matrix with less than 5 columns');
 return
end
result=zeros(m,n);
result(:,3:n-2)=(arr(:,5:n)-8*arr(:,4:n-1)+8*arr(:,2:n-3)-arr(:,1:n-4)) ...
               ./(yy(:,4:n-1)-yy(:,2:n-3))/6;
result(:,1)=result(:,3);
result(:,2)=result(:,3);
result(:,n-1)=result(:,n-2);
result(:,n)=result(:,n-2);
