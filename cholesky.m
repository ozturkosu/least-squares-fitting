% 2019862s 
% Assessed Coursework 1: Least Squares Fitting

function [R,A] = cholesky(A)

% Associate A with R
R=A;

% Store dimensions of A
[m,n]=size(A);

% Loop over rows
for k=1:m,
    % Nested loop over rows
    for j=(k+1):m
        % Formula from lecture notes
        R(j,(j:m))=R(j,(j:m))-(R(k,(j:m))*(conj(R(k,j))))...
            /(R(k,k));
    end
    R(k,(k:m))=R(k,(k:m))/(sqrt(R(k,k)));
    R(k,(1:k-1))=0;
end
% Extract reduced R
R=R(1:n,:);
end

