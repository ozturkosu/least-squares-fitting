% 2019862s 
% Assessed Coursework 1: Least Squares Fitting

%% Reduced QR via Householder triangularization
% Provide matrix A in the console
A=input('Enter the Matrix A: \n');
[houseQ,houseR] = reducedQRHouseholder(A);



%% Least Squares Question 5 - Run this to compare three methods
format long 
m = 50;
n = 12;
t = linspace(0,1,m);

V = vander(t);
V = fliplr(V);
A = zeros(m,n);
    for i=1:n,
        A(:,i) = V(:,i); 
    end
b = sin(5*t')

%%%%%%%%%%%%% QR via Householder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute reduced QR factorization by Householder
[houseQ, houseR] = reducedQRHouseholder(A);
% Conjugate transpose of reducedQ
conjQ = (conj(houseQ)');

% Perform backward substitution
householdery = (conjQ)*b;
householderx = backSubstitution(houseR,householdery)

%%%%%%%%%%%%% Gram-Schmidt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute QR factorization by Gram-Schmidt
[gramSchmidtQ,gramSchmidtR]=modifiedGramSchmidt(A);

% Conjugate transpose of reducedQ
conjGramSchmidtQ = (conj(gramSchmidtQ'));

% Perform backward substitution
gramSchmidty = (conjGramSchmidtQ)*b;
gramSchmidtx = backSubstitution(gramSchmidtR,gramSchmidty)

%%%%%%%%%%%%% Normal Equations via Cholesky %%%%%%%%%%%%%%%%%%%%%%
% Compute Aconj(A) and conj(A)b
prodA = (conj(A'))*A;
choleskyy = (conj(A'))*b;

% Compute Cholesky factorization of A
[choleskyR] = (cholesky(prodA));

% Solve lower triangular system using forward substitution
w = forwSubstitution((conj(choleskyR')),choleskyy);

choleskyx = backSubstitution(choleskyR,w)

%% Absolute Error
m = 50;
t = linspace(0,1,m);
absoluteErrorHouseholder = abs(A*householderx-b);
absoluteErrorGramSchmidt = abs(A*gramSchmidtx-b);
absoluteErrorCholesky = abs(A*choleskyx-b);

fig1=figure
    plot(t, absoluteErrorHouseholder,'r');
    hold on
    plot(t, absoluteErrorGramSchmidt,'k--');
    plot(t, absoluteErrorCholesky,'k');
    title('Absolute Errors vs. t for m=50, n=12');
    xlabel('t');
    ylabel('Absolute Error for each Method');
    legend('Householder','Gram-Schmidt','Cholesky','Location','northwest');
    print(fig1,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/absoluteError.jpeg','-djpeg');
    figure(fig1);
    
fig2=figure
    semilogy(t, absoluteErrorHouseholder,'r');
    hold on
    semilogy(t, absoluteErrorGramSchmidt,'k--');
    semilogy(t, absoluteErrorCholesky,'k');
    title('Absolute Errors in log scale vs. t');
    xlabel('t');
    ylabel('Absolute Error in log scale for each Method');
    legend('Householder','Gram-Schmidt','Cholesky','Location','northwest');
    print(fig2,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/absoluteErrorLog.jpeg','-djpeg');    
    figure(fig2);

%% Relative Error

format long 
m = 50;
t = linspace(0,1,m);
b = sin(5*t');
V = vander(t);
V = fliplr(V);
for n = 5:30,
        A = V(:,1:n);
        

%%%%%%%%%%%%% QR via Householder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute reduced QR factorization by Householder
[houseQ, houseR] = reducedQRHouseholder(A);
% Conjugate transpose of reducedQ
conjQ = (conj(houseQ)');

% Perform backward substitution
householdery = (conjQ)*b;
householderx = backSubstitution(houseR,householdery);

%%%%%%%%%%%%% QR via Gram-Schmidt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute QR factorization by Gram-Schmidt
[gramSchmidtQ,gramSchmidtR]=modifiedGramSchmidt(A);

% Conjugate transpose of reducedQ
conjGramSchmidtQ = (conj(gramSchmidtQ'));

% Perform backward substitution
gramSchmidty = (conjGramSchmidtQ)*b;
gramSchmidtx = backSubstitution(gramSchmidtR,gramSchmidty);

%%%%%%%%%%%%% Normal Equations via Cholesky %%%%%%%%%%%%%%%%%%%%%%
% Compute Aconj(A) and conj(A)b
prodA = (conj(A'))*A;
choleskyy = (conj(A'))*b;

% Compute Cholesky factorization of A
[choleskyR] = (cholesky(prodA));

% Solve lower triangular system using forward substitution
w = forwSubstitution((conj(choleskyR')),choleskyy);

choleskyx = backSubstitution(choleskyR,w);

    relativeErrorHouseholder(n-4) = (norm((A*householderx-b),2))/(norm(b,2))
    relativeErrorGramSchmidt(n-4) = (norm((A*gramSchmidtx-b),2))/(norm(b,2))
    relativeErrorCholesky(n-4) = (norm((A*choleskyx-b),2))/(norm(b,2))
end
fig3=figure
    plot(linspace(5,30,26), relativeErrorHouseholder,'r');
    hold on
    plot(linspace(5,30,26), relativeErrorGramSchmidt,'k--');
    plot(linspace(5,30,26), relativeErrorCholesky,'k');
    title('Relative Errors vs. n=[5,30]');
    xlabel('n=[5,30]');
    ylabel('Relative Error for each Method');
    legend('Householder','Gram-Schmidt','Cholesky','Location','northwest');
    print(fig3,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/relativeError.jpeg','-djpeg');

fig4=figure
    semilogy(linspace(5,30,26), relativeErrorHouseholder,'r');
    hold on
    semilogy(linspace(5,30,26), relativeErrorGramSchmidt,'k--');
    semilogy(linspace(5,30,26), relativeErrorCholesky,'k');
    title('Relative Errors in log scale vs. n=[5,30].');
    xlabel('n=[5,30]');
    ylabel('Relative Error in log scale for each Method');
    legend('Householder','Gram-Schmidt','Cholesky','Location','northwest');
    print(fig4,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/relativeErrorLog.jpeg','-djpeg');

%% Condition number
format long 
m = 50;
t = linspace(0,1,m);
b = sin(5*t');
V = vander(t);
V = fliplr(V);
for n = 5:30,
    A = V(:,1:n);
    conditionNumberA(n-4)=cond(A);
    conditionNumberAStarA(n-4)=cond((conj(A))'*A);
end

fig5=figure
    plot(linspace(5,30,26), conditionNumberA,'r');
    hold on
    plot(linspace(5,30,26), conditionNumberAStarA,'k--');
    title('Condition number of A vs. n=[5,30].');
    xlabel('n=[5,30]');
    ylabel('Condition Number of A for each Method');
    legend('Condition Number of A','Condition Number of A*A','Location','northwest');
    print(fig5,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/conditionNumber.jpeg','-djpeg');

fig6=figure
    semilogy(linspace(5,30,26), conditionNumberA,'r');
    hold on
    semilogy(linspace(5,30,26), conditionNumberAStarA,'k--');
    title('Condition number of A in log scale vs. n=[5,30].');
    xlabel('n=[5,30]');
    ylabel('Condition Number of A for each Method');
    legend('Condition Number of A','Condition Number of A*A','Location','northwest');
    print(fig6,'/Users/yanastaneva/Library/Mobile Documents/com~apple~CloudDocs/Year 5/Numerical Methods/Assignments/AC1/conditionNumberLog.jpeg','-djpeg');
