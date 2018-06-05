% 2019862s 
% Assessed Coursework 1: Least Squares Fitting

function x=backSubstitution(U,b)

% Store the dimensions of the upper triangular U
[m,n]=size(U);

% Initiate the zero colum vector 
x=zeros(m,1);

% Iterate over the rows
for j=m:-1:1
   % Compute the j-th entry of x
   x(j) = ( b(j) - U(j, :)*x )/U(j, j);
end