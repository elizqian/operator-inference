function x = tikhonov(b,A,k)
% solves linear regression Ax = b with L2/tikhonov regularization penalty
%
% INPUTS
% A     data matrix
% b     right-hand side
% k     Tikhonov weighting
%
% OUTPUT
% x     solution
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

[~,q] = size(b);
[~,p] = size(A);
pseudo = sqrt(k) * eye(p);
Aplus  = [A;pseudo];
bplus  = [b;zeros(p,q)];

x = Aplus\bplus;
end