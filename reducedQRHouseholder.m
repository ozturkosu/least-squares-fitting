% 2019862s 
% Assessed Coursework 1: Least Squares Fitting

function [Q,R] = reducedQRHouseholder(A)

% Store the dimensions of A.
[m, n] = size(A);

% Associate the matrix R with the input A initially.
R = A; 
% Initialize an mxn unitary matrix Q, all entries are zeros.
Q = eye(m,n);
% Initialize matrix W to store vectors v
V = zeros(m,n);
% For loop over the columns.
for j = 1:n,
    
  % Take j-th column of R.
  x = R(j:m,j);
  
  % Set sign(x(1))=1 if x(1)=0.
  if x(1)==0, 
     s=1;
  % Otherwise use formula as usual.
  else 
     s=sign(x(1));
  end
  
  % Form the reflection vector.
  x(1)=s*norm(x,2)+x(1);
  
  % Associate vector v with x.
  v=x;
  
  % Normalize the reflection vector.
  v = v./norm(v,2);
  
  % Store the vectors in a matrix to be later used for
  % implicit calculation of reducedQ.
  V(j:m,j)=v;
  
  % Form the upper triangular R by clearing entries below
  % the diagonal in the j-th column at the j-th iteration.
  % This is equivalent to left multiplying by the Q_k matrix.
  R(j:m, j:n) = R(j:m, j:n) - ...
      2*v*(v'*R(j:m, j:n));
end
  % Implicitly calculate the product reducedQ*v by extracting 
  % v vectors from the matrix V
  for i=1:n,
      for k=n:-1:1;
        Q(k:m,i) = Q(k:m,i) - 2*V(k:m,k)*...
            (V(k:m,k)'*Q(k:m,i));
      end
  end
  % Ensure to extract only non-zero rows from R
  R=R(1:n,:);
end