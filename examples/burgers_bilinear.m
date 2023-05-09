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

clear; clc;
addpath('../')

%% Problem set-up
n       = 2^7;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mu = 0.1;               % diffusion coefficient

% run FOM with input 1s to get reference trajectory
u_ref = ones(K,1);
% u_ref = sin(pi/2 * linspace(0,1,K))';

[A,B,F,N] = getBurgersMatrices(n,1/(n-1),mu);
s_ref = semiImplicitEuler(A,F,B,N,dt,u_ref);

%% Surface Plot for verification
if true
    s = surf(linspace(0.0,T_end,K+1),linspace(0.0,1.0,n),s_ref,'FaceAlpha',0.8);
    s.EdgeColor = 'none';
    xlabel("t");
    ylabel("x");
end

%% Operator inference parameters
params.modelform = 'LQIB';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

%% collect data for a series of trajectories with random inputs
num_inputs = 10;
U_rand = rand(K,num_inputs);
x_all = cell(num_inputs,1);
xdot_all = cell(num_inputs,1);
for i = 1:num_inputs
    s_rand = semiImplicitEuler(A,F,B,N,dt,U_rand(:,i));
    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:30;
err_inf = zeros(length(r_vals),1);
err_int = zeros(length(r_vals),1);

% intrusive
Vr = U_svd(:,1:max(r_vals));
Aint = Vr' * A * Vr;
Bint = Vr' * B;
Ln = elimat(n); Dr = dupmat(max(r_vals));
Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;
Nint = Vr' * N * Vr;

% op-inf
[operators] = inferOperators(X, U, Vr, params, R);
Ahat = operators.A;
Fhat = operators.F;
Bhat = operators.B;
Nhat = operators.N;

for j = 1:length(r_vals)
    r = r_vals(j);
    Vr = U_svd(:,1:r);

    Fhat_extract = extractF(Fhat, r);
    s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),Nhat(1:r,1:r),dt,u_ref);
    s_rec = Vr*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
    
    Fint_extract = extractF(Fint, r);
    s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),Nint(1:r,1:r),dt,u_ref);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
end


%% Plotting
figure(2); clf;
semilogy(r_vals,err_inf, DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals,err_int, DisplayName="int"); 
hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error','Interpreter','LaTeX')


%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler(Ahat, Fhat, Bhat, Nhat, dt, u_input)
    K = size(u_input,1);
    r = size(Ahat,1);
    s_hat = zeros(r,K+1); % initial state is zeros everywhere
    ImdtA = eye(r) - dt*Ahat;

    for i = 1:K
        ssq = get_x_sq(s_hat(:,i)')';
        s_hat(:,i+1) = ImdtA\(s_hat(:,i) + dt*Fhat*ssq + dt*Bhat*u_input(i) + dt*(Nhat*s_hat(:,i))*u_input(i));
        if any(isnan(s_hat(:,i)))
            warning(['ROM unstable at ',num2str(i-1),'th timestep'])
            break
        end
    end
end

%% builds matrices for Burgers full-order model
function [A, B, F, N] = getBurgersMatrices(n,dx,mu)
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