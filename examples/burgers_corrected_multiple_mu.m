% Elizabeth Qian (elizqian@mit.edu) 11 July 2019
% Tomoki Koike (tkoike3@gatech.edu) 10 May 2023 [EDITED]
%
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
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mus = 0.1:0.1:1.0; % diffusion coefficient
M = length(mus);

% Operator inference parameters
params.modelform = 'LQI';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

% run FOM with input 1s to get reference trajectory
u_ref = ones(K,1);

% Settings
r_vals = 1:18;
err_inf = zeros(length(mus),length(r_vals));
err_int = zeros(length(mus),length(r_vals));

% Store all s_ref values for verification
s_ref_all = cell(M,1);

% Store all operators to use for test data
infop_all = cell(M,1);
intop_all = cell(M,1);

%% LEARN AND ANALYZE TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:M
    mu = mus(i);
    [A,B,F] = getBurgersMatrices(N,1/(N-1),dt,mu);
    s_ref = semiImplicitEuler(A,F,B,dt,u_ref);
    s_ref_all{i} = s_ref;

    %% collect data for a series of trajectories with random inputs
    num_inputs = 10;
    U_rand = rand(K,num_inputs);
    x_all = cell(num_inputs,1);
    xdot_all = cell(num_inputs,1);
    for k = 1:num_inputs
        s_rand = semiImplicitEuler(A,F,B,dt,U_rand(:,k));
        x_all{k}    = s_rand(:,2:end);
        xdot_all{k} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
    end
    
    X = cat(2,x_all{:});  % concatenate data from random trajectories
    R = cat(2,xdot_all{:});    
    U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);
    [U_svd,~,~] = svd(X,'econ');  % take SVD for POD basis

    % intrusive
    Vr = U_svd(:,1:max(r_vals));
    Aint = Vr' * A * Vr;
    Bint = Vr' * B;
    Ln = elimat(N); Dr = dupmat(max(r_vals));
    Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;
    op_int.A = Aint; op_int.B = Bint; op_int.F = Fint;
    intop_all{i} = op_int;
    
    % op-inf
    [operators] = inferOperators(X, U, Vr, params, R);
    Ahat = operators.A;
    Fhat = operators.F;
    Bhat = operators.B;
    infop_all{i} = operators;

    %% For different basis sizes r, compute basis, learn model, and calculate state error 
    for j = 1:length(r_vals)
        r = r_vals(j);
        Vr = U_svd(:,1:r);
        
        % Extract operators for inferred model
        Fhat_extract = extractF(Fhat, r);
        s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref);
        s_rec = Vr*s_hat;
        err_inf(i,j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
        
        % Extract operators for intrusive model
        Fint_extract = extractF(Fint, r);
        s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref);
        s_tmp = Vr*s_int;
        err_int(i,j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
    end
end

%% Plot relative state error
err_inf_avg = mean(err_inf(2:end,:));
err_int_avg = mean(err_int(2:end,:));

figure(1); clf;
semilogy(r_vals(1:14),err_inf_avg(1:14), DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals(1:14),err_int_avg(1:14), DisplayName="int"); hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error (Train)','Interpreter','LaTeX')

%% Verify FOM
figure(2)
t = tiledlayout(2,5);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for i = 1:M
    nexttile
        s = surf(linspace(0.0,T_end,K+1),linspace(0.0,1.0,N),s_ref_all{i},'FaceAlpha',0.8,DisplayName="$\mu$="+mus(i));
        s.EdgeColor = 'none';
        xlabel("t");
        ylabel("x");
        legend(Interpreter="latex");
end


%% ANALYZE TEST DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mus_new = rand(5,1);
M_new = length(mus_new);

err_inf_new = zeros(M_new,length(r_vals));
err_int_new = zeros(M_new,length(r_vals));

%%
for i = 1:M_new
    mu = mus_new(i);
    [A,B,F] = getBurgersMatrices(N,1/(N-1),dt,mu);
    s_ref = semiImplicitEuler(A,F,B,dt,u_ref);

    % Interpolate the intrusive and inferred operator for random mu value
    [Aint, Bint, Fint] = spline_operators(intop_all, max(r_vals), mus, mu);
    [Ainf, Binf, Finf] = spline_operators(infop_all, max(r_vals), mus, mu);

    %% For different basis sizes r, compute basis, learn model, and calculate state error 
    for j = 1:length(r_vals)
        r = r_vals(j);
        Vr = U_svd(:,1:r);
        
        % Extract operators for inferred model
        Fhat_extract = extractF(Fhat, r);
        s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref);
        s_rec = Vr*s_hat;
        err_inf_new(i,j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
        
        % Extract operators for intrusive model
        Fint_extract = extractF(Fint, r);
        s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref);
        s_tmp = Vr*s_int;
        err_int_new(i,j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
    end
end

%% Plotting
err_inf_avg_new = mean(err_inf_new);
err_int_avg_new = mean(err_int_new);

figure(3); clf;
semilogy(r_vals,err_inf_avg_new, DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals,err_int_avg_new, DisplayName="int"); hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error (Test)','Interpreter','LaTeX')


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

%% builds matrices for Burgers full-order model
function [A, B, F] = getBurgersMatrices(N,dx,dt,mu)
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