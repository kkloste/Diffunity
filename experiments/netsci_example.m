load netscience-cc.mat;

% Get eigenmatrix
d = full(sum(A,1));
dsqinv = diag( d.^(-1/2) );
Anormalized = dsqinv*A*dsqinv;
Anormalized = full(Anormalized);

[V,Lambda] = eig(Anormalized);
lambdas = diag(Lambda);
L = diag(d) - A;

n = size(A,1);


%%

s = randi(n); % get seed.

ejd = sparse( s, 1, d(s)^0.5 , n, 1);
temp = V'*ejd;
h = temp.*temp;
g = h - Lambda*h;

cvx_begin quiet
    cvx_precision high
    variable x(n);
    minimize quad_over_lin(g'*x,h'*x);
    subject to
        x >= 0;   % check y >= 0
        %x <= 1; % check y(S) >= 1
cvx_end


%%

P = polyfit(lambdas,x,20);


format long g
pred = polyval(P,lambdas) - x;

diffP = hpolym( A*dsqinv*dsqinv , P);
fdiff = diffP(:,s);


