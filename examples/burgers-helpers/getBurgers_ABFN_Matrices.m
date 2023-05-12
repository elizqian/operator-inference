%% builds matrices for Burgers full-order model
function [A, B, F, N] = getBurgers_ABFN_Matrices(n,dx,mu)
    % construct linear operator resulting from second derivative
    A = mu*gallery('tridiag',n,1,-2,1)/(dx^2);

    % construct quadratic operator F using central difference for the derivative
    ii = reshape(repmat(2:n-1,2,1),2*n-4,1);
    m = 2:n-1;
    mi = n*(n+1)/2 - (n-m).*(n-m+1)/2 - (n-m);              % this is where the xi^2 term is
    mm = n*(n+1)/2 - (n-m).*(n-m+1)/2 - (n-m) - (n-(m-2));  % this is where the x_{i-1}^2 term is
    jp = mi + 1;        % this is the index of the x_{i+1}*x_i term
    jm = mm + 1;        % this is the index of the x_{i-1}*x_i term
    jj = reshape([jp; jm],2*n-4,1);
    vv = reshape([ones(1,n-2); -ones(1,n-2)],2*n-4,1)/(2*dx);
    F = -sparse(ii,jj,vv,n,n*(n+1)/2);
    F(1,2) = -1/2/dx; F(n,end-1) = 1/2/dx;

    % construct input matrix B
    B = [1; zeros(n-2,1); -1] * mu / dx^2;

    % construct bilinear matrix N
    N = zeros(n,n);
    N(1,1) = 1/2/dx;
    N(n,n) = 1/2/dx;  % it's supposed to be (-1) but since the input is -u(t) at the boundary, we multiply this by (-1)
end