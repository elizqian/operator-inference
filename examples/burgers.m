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
[Uref_svd,~,~] = svd(s_ref,"econ");

%% Operator inference parameters
params.modelform = 'LQI';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

%% collect data for a series of trajectories with random inputs
num_inputs = 10;
U_rand = rand(K,num_inputs);
x_all = cell(num_inputs,1);
xdot_all = cell(num_inputs,1);
for i = 1:num_inputs
    s_rand = burgersFOM(N,dt,T_end,mu,U_rand(:,i));
    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%%
tdata = 0:dt:T_end;
xdata = linspace(0,1,N);
surf(tdata(2:end),xdata,x_all{1},EdgeColor="none",FaceAlpha=0.8);

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:10;
err_inf = zeros(length(r_vals),1);
err_int = zeros(length(r_vals),1);
for j = 1:length(r_vals)
    r = r_vals(j);
    Vr = U_svd(:,1:r);
    
    [operators] = inferOperators(X, U, Vr, params, R);
    Ahat = operators.A;
    Fhat = operators.F;
    Bhat = operators.B;

    s_hat = semiImplicitEuler(Ahat,Fhat,Bhat,dt,u_ref);
    s_rec = Vr(:,1:r)*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')^2/norm(s_ref,'fro')^2;

    % intrusive
    rr = r;
    Ur = Uref_svd(:,1:rr);
    Aint = Ur' * A * Ur;
    Bint = Ur' * B;
    Ln = elimat(N); Dr = dupmat(rr);
    Fint = Ur' * F * Ln * kron(Ur,Ur) * Dr;

    s_int = semiImplicitEuler(Aint,Fint,Bint,dt,u_ref);
    s_tmp = Ur(:,1:rr)*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')^2/norm(s_ref,'fro')^2;
end

figure(2); clf
semilogy(r_vals,err_inf); grid on
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error','Interpreter','LaTeX')

figure(3); clf
semilogy(r_vals,err_int); grid on
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers intrusive model error','Interpreter','LaTeX')

%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler(Ahat, Fhat, Bhat, dt, u_input)
K = size(u_input,1);
r = size(Ahat,1);
s_hat = zeros(r,K+1); % initial state is zeros everywhere
ImdtA = eye(r) - dt*Ahat;
for i = 1:K
    ssq = get_x_sq(s_hat(:,i)')';
    s_hat(:,i+1) = ImdtA\(s_hat(:,i) + dt*Fhat*ssq + dt*Bhat*u_input(i));
    if any(isnan(s_hat(:,i+1)))
        warning(['ROM unstable at ',num2str(i),'th timestep'])
        break
    end
end
end

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

%% Other
function D2 = dupmat(n)
  m   = n * (n + 1) / 2;
  nsq = n^2;
  r   = 1;
  a   = 1;
  v   = zeros(1, nsq);
  cn  = cumsum(n:-1:2);   % [EDITED, 2021-08-04], 10% faster
  for i = 1:n
     % v(r:r + i - 2) = i - n + cumsum(n - (0:i-2));
     v(r:r + i - 2) = i - n + cn(1:i - 1);   % [EDITED, 2021-08-04]
     r = r + i - 1;
     
     v(r:r + n - i) = a:a + n - i;
     r = r + n - i + 1;
     a = a + n - i + 1;
  end
  
  D2 = sparse(1:nsq, v, 1, nsq, m);
end

function L = elimat(m)
  T = tril(ones(m)); % Lower triangle of 1's
  f = find(T(:)); % Get linear indexes of 1's
  k = m*(m+1)/2; % Row size of L
  m2 = m*m; % Colunm size of L
  L = zeros(m2,k); % Start with L'
  x = f + m2*(0:k-1)'; % Linear indexes of the 1's within L'
  L(x) = 1; % Put the 1's in place
  L = L'; % Now transpose to actual L
end