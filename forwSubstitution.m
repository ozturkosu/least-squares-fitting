% 2019862s 
% Assessed Coursework 1: Least Squares Fitting

function w=forwSubstitution(L,b)

% Store the dimensions of the lower triangular L
[m,n]=size(L);

% Initiate a column vector 
n = length(b);

% Initiate a matrix for the output
w = zeros(m,1);
    % Loop over rows
    for i=1:m,
        % Compute the i-th entry from the substitution
        w(i) = ( b(i) - L(i, :)*w )/L(i, i);
    end

end 


