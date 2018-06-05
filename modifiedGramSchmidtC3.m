function [Q,R]=modifiedGramSchmidtC3(A)

% Store dimensions of matrix A
[m,n]=size(A);

% Inititate zero matrices for the reduced Q and R
Q=zeros(m,n);
R=zeros(m,n);

% Initiate a zero column vector 
v=zeros(m,1);

% For loop for extracting columns from the matrix A
for i=1:n;
    v(:,i)=A(:,i);
end

% For loop which implements the modified Gram-Schmidt
for i=1:n;
    % Diagonal entry of R is the 2-norm of the i-th column vector 
    R(i,i)=norm(v(:,i),2);
    % Compute the i-th column of Q 
    Q(:,i)=v(:,i)/R(i,i);
       for j=i+1:n;
           % Compute the off-diagonal entries of R
           R(i,j)=Q(:,i)'*v(:,j);
           % Update the next column vector to be used 
           % in the next iteration
           v(:,j)=v(:,j)-R(i,j)*Q(:,i);
       end
end

% Make sure to extract the reduced R
R=R(1:n,:);

        
    
    