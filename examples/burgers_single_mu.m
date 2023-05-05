%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
addpath('../')

%% Problem set-up
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mu = 0.1;

% run FOM with input 1s to get reference trajectory
u_ref = ones(K,1);

[A,B,F] = getBurgersMatrices(N,1/(N-1),mu);
s_ref = semiImplicitEuler(A,F,B,dt,u_ref);

%% Surface Plot for verification
if true
    s = surf(linspace(0.0,T_end,K),linspace(0.0,1.0,N),s_ref,'FaceAlpha',0.8);
    s.EdgeColor = 'none';
    xlabel("t");
    ylabel("x");
end

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
    s_rand = semiImplicitEuler(A,F,B,dt,U_rand(:,i));
    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end
X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);
[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:15;
err_inf = zeros(length(r_vals),1);
err_int = zeros(length(r_vals),1);
diff = zeros(length(r_vals),1);

% intrusive
Vr = U_svd(:,1:max(r_vals));
Aint = Vr' * A * Vr;
Bint = Vr' * B;
Ln = elimat(N); Dr = dupmat(max(r_vals));
Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;

% op-inf
[operators] = inferOperators(X, U, Vr, params, R);
Ahat = operators.A;
Fhat = operators.F;
Bhat = operators.B;

for j = 1:length(r_vals)
    r = r_vals(j);
    Vr = U_svd(:,1:r);
    
    Fhat_extract = extractF(Fhat, r);
    s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref);
    s_rec = Vr*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
    
    Fint_extract = extractF(Fint, r);
    s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
end

figure(2); clf
semilogy(r_vals,err_inf, DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals,err_int, DisplayName="int"); 
hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error','Interpreter','LaTeX')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler(Ahat, Fhat, Bhat, dt, u_input)
    K = size(u_input,1);
    r = size(Ahat,1);
    s_hat = zeros(r,K); % initial state is zeros everywhere
    ImdtA = eye(r) - dt*Ahat;

    for i = 2:K
        ssq = get_x_sq(s_hat(:,i-1)')';
        if i == 2
            Du = u_input(i) - 0.0;
        else
            Du = u_input(i) - u_input(i-1);
        end
        s_hat(:,i) = ImdtA\(s_hat(:,i-1) + dt*Fhat*ssq + Bhat*Du);
        if any(isnan(s_hat(:,i)))
            warning(['ROM unstable at ',num2str(i-1),'th timestep'])
            break
        end
    end
end

%% solves Burgers equation from zero initial condition with specified input
% function [s_all,A,B,F] = burgersFOM(N,dt,T_end,mu,u)
%     dx = 1/(N-1);
%     
%     K = T_end/dt;
%     
%     [A,B,F] = getBurgersMatrices(N,dx,mu);
%     ImdtA = eye(N)-dt*A;
%     ImdtA(1,1:2) = [1 0]; ImdtA(N,N-1:N) = [0 1]; % Dirichlet boundary conditions
%     
%     s_all = zeros(N,K);       % initial state is zero everywhere
%     for i = 2:K
%         ssq = get_x_sq(s_all(:,i-1)')';
%         s_all(:,i) = ImdtA\(s_all(:,i-1) + dt*F*ssq + dt*B*u(i-1));
%     end
% end

%% builds matrices for Burgers full-order model
function [A, B, F] = getBurgersMatrices(N,dx,mu)
    % construct linear operator resulting from second derivative
    A = mu*gallery('tridiag',N,1,-2,1)/(dx^2);
    A(1,1:2) = [0 0]; A(N,N-1:N) = [0 0]; % Dirichlet boundary conditions

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
    F = -sparse(ii,jj,vv,N,N*(N+1)/2);  % CHANGE: MULTIPLIED BY (-1)

    % construct input matrix B
    B = [1; zeros(N-2,1); -1];
end
