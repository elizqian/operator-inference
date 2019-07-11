% Sets up operator inference problem for 1D heat equation as described in
%   Peherstorfer, B. and Willcox, K., "Data-driven operator inference for
%   non-intrusive projection-based model reduction." Computer Methods in
%   Applied Mechanics and Engineering, 306:196-215, 2016.

clear
addpath('../')

%% Problem set-up
T       = 1;      % final time
dx      = 2^-7;   % mesh width
dt      = 1e-3;   % time step
x_grid  = 0:dx:1; % grid points
t       = 0:dt:T; % time points

N = length(x_grid)-2;

U = ones(length(t),1);  % Dirichlet BCs (symmetric)
x0 = zeros(N,1);            % initial condition

M = 10;                     % number of diffusivity values between [0.1, 10] to learn models for
mu_all = logspace(-1,1,M);  % diffusivity values

%% Operator inference parameters
params.modelform = 'LI';            % model is linear with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = 'BE';            % scheme for estimating time derivative

%% Infer operators for reduced models for each value of diffusivity
all_Ahat = cell(M,1);
all_Bhat = cell(M,1);
all_X = cell(M,1);
all_Vr = cell(M,1);

r_max = 8;          % maximum model size
for i = 1:M
    % run FOM for specified mu to collect data
    mu = mu_all(i);
    [A,B] = getHeatMatrices(N,dx,mu);
    X = backwardEuler(x0,A,B,U,t);
    
    % compute POD basis for FOM data
    [Usvd,~,~] = svd(X);
    Vr = Usvd(:,1:r_max);
    all_X{i} = X;
    all_Vr{i} = Vr;
    
    % infer operators from data
    operators = inferOperators(X, U, Vr, params);
    all_Ahat{i} = operators.A;
    all_Bhat{i} = operators.B;
end

%% calculate state and projection errors for different basis sizes
proj_errors = zeros(r_max,1);
state_errors = zeros(r_max,1);
for n = 1:r_max         % loop through basis sizes
    for i = 1:M     % loop through different diffusivities
        X  = all_X{i};
        Vr = all_Vr{i}(:,1:n);
        proj_errors(n) = proj_errors(n) + 1/M*norm(X - Vr*Vr'*X,'fro')^1/norm(X,'fro')^1;
        
        Ahat = all_Ahat{i};
        Bhat = all_Bhat{i};
        XNhat = backwardEuler(Vr'*x0,Ahat(1:n,1:n),Bhat(1:n,:),U,t);
        state_errors(n) = state_errors(n) + 1/M*norm(X-Vr*XNhat,'fro')^1/norm(X,'fro')^1;
    end
end

%% plot errors
% plot projection errors at different basis sizes
figure(1); clf
semilogy(proj_errors,'b','LineWidth',2); grid on; hold on
semilogy(proj_errors,'bo','MarkerSize',10,'LineWidth',2);
xlabel('dimension n')
ylabel('Projection error (eq. 29 w/o square)')
title('Projection error, heat equation, training set')

% plot state errors at different basis sizes
figure(2); clf
semilogy(state_errors,'b','LineWidth',2); grid on; hold on
semilogy(state_errors,'bo','MarkerSize',10,'LineWidth',2);
xlabel('dimension n')
ylabel('State error (eq. 30 w/o square)')
title('State error, heat equation, training set')

%% Function to build FOM matrices
function [A,B] = getHeatMatrices(N,dx,mu)
    A = mu*(-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/dx^2;
    B = [mu/dx^2; zeros(N-2,1); mu/dx^2];
end

%% Backward Euler integrator
function state = backwardEuler(x0,A,B,u_all,time)
K = length(time);
N = length(x0);
state = zeros(N,K);
state(:,1) = x0;
for i = 2:K 
    u = u_all(i-1);
    dt = time(i) - time(i-1);
    nextstate = (A - eye(N)/dt)\(-state(:,i-1)/dt - B*u);
    state(:,i) = nextstate;
end
end