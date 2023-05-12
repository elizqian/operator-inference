%% builds matrices for Burgers full-order model
function [A, B, F] = getBurgers_ABF_Matrices(N,dx,dt,mu)
    % construct linear operator resulting from second derivative
    A = mu*gallery('tridiag',N,1,-2,1)/(dx^2);
    A(1,1:2) = [-1/dt 0]; A(N,N-1:N) = [0 -1/dt]; % Dirichlet boundary conditions

    % construct quadratic operator F using central difference for the
    % derivative
    ii = reshape(repmat(2:N-1,2,1),2*N-4,1);
    m = 2:N-1;
    mi = N*(N+1)/2 - (N-m).*(N-m+1)/2 - (N-m);              % this is where the xi^2 term is
    mm = N*(N+1)/2 - (N-m).*(N-m+1)/2 - (N-m) - (N-(m-2));  % this is where the x_{i-1}^2 term is
    jp = mi + 1;        % this is the index of the x_{i+1}*x_i term
    jm = mm + 1;        % this is the index of the x_{i-1}*x_i term
    jj = reshape([jp; jm],2*N-4,1);
    vv = reshape([ones(1,N-2); -ones(1,N-2)],2*N-4,1)/(2*dx);
    F = -sparse(ii,jj,vv,N,N*(N+1)/2);

    % construct input matrix B
    B = [1; zeros(N-2,1); -1]/dt;
end