% Elizabeth Qian (elizqian@mit.edu) 11 July 2019
% -----------------------------------------------------------------------
% Based on operator inference problem for Burgers equation described in:
%   Peherstorfer, B. and Willcox, K., "Data-driven operator inference for
%   non-intrusive projection-based model reduction." Computer Methods in
%   Applied Mechanics and Engineering, 306:196-215, 2016.
%
% Note this script only considers one value of the parameter mu (whereas
% the above paper considers multiple values)
%
% See also:
%   Qian, E., Kramer, B., Marques, A. and Willcox, K., "Transform & Learn:
%   A data-driven approach to nonlinear model reduction." In AIAA Aviation 
%   2019 Forum, June 17-21, Dallas, TX.

clear
addpath('../')

%% Problem set-up
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mu = 0.1;               % diffusion coefficient

% run FOM with input 1s to get reference trajectory
u_ref = ones(K,1);
[s_ref,A,B,F] = burgersFOM(N,dt,T_end,mu,u_ref);
H = F2H(F);

%%
% check index-wise constraint
get_h = @(i,j,k) H(i,(k-1)*N+j);
derp = zeros(N,N,N);
derpsum = 0;
for i = 1:N
    for j = 1:N
        for k = 1:N
            foo = get_h(i,j,k) + get_h(j,i,k) + get_h(k,j,i);
            derp(i,j,k) = get_h(i,j,k) + get_h(j,i,k) + get_h(k,j,i);
            derpsum = derpsum + abs(foo);
        end
    end
end

% check inner product
viol = 0;
for i = 1:10
    randtest = rand(N,1);
    viol = viol + abs(randtest'*H*kron(randtest,randtest));
end
viol

%% solves Burgers equation from zero initial condition with specified input
function [s_all,A,B,F] = burgersFOM(N,dt,T_end,mu,u)
dx = 1/(N-1);

K = T_end/dt;

[A,B,F] = getBurgersMatrices(N,dx,mu);
ImdtA = eye(N)-dt*A;
ImdtA(1,1:2) = [1 0]; ImdtA(N,N-1:N) = [0 1]; % Dirichlet boundary conditions

s_all = zeros(N,K+1);       % initial state is zero everywhere
for i = 1:K
    ssq = get_x_sq(s_all(:,i)')';
    s_all(:,i+1) = ImdtA\([0; s_all(2:N-1,i); 0] - dt*F*ssq + dt*B*u(i));
end
end

%% builds matrices for Burgers full-order model
function [A, B, F] = getBurgersMatrices(N,dx,mu)
    % construct linear operator resulting from second derivative
    A = mu*gallery('tridiag',N,1,-2,1)/(dx^2);
    A(1,1:2) = [1 0]; A(N,N-1:N) = [0 1]; % Dirichlet boundary conditions

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
    F = sparse(ii,jj,vv,N,N*(N+1)/2);

    % construct input matrix B
    B = [1; zeros(N-2,1); -1];
end