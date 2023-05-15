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

clear; close all; clc;
addpath('../',"burgers-helpers/");

%% Problem set-up
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mu = 0.5;               % diffusion coefficient

type = 3;
% run FOM with input 1s to get reference trajectory
if type == 1
    u_ref = ones(K,1);
    IC = zeros(N,1);
elseif type == 2
    u_ref = zeros(K,1);
    IC = -exp(linspace(0,1,N) - 1)' .* sin(pi * linspace(0,1,N) - pi)';
else
    u_ref = zeros(K,1);
    IC = exp(-10*(2*linspace(0,1,N) - 1).^2)';
    fic = @(a) exp(-10*(a*linspace(0,1,N) - 1).^2);
    ic_a = linspace(3,7.5,10);
end

[A, B, F] = getBurgers_ABF_Matrices(N,1/(N-1),dt,mu);
H = F2Hs(F);
s_ref = semiImplicitEuler(A, F, B, dt, u_ref, IC);

%% Surface Plot for verification
figure(1);
s = surf(linspace(0.0,T_end,K+1),linspace(0.0,1.0,N),s_ref,'FaceAlpha',0.8,DisplayName="\mu="+num2str(mu));
s.EdgeColor = 'none';
xlabel("t, time");
ylabel("\omega, space");
zlabel("x(\omega,t), velocity")
axis tight
view(-73.25,38.649)
grid on
text(0.1,0.8,1,"\mu = "+num2str(mu),'FontSize',14);

%% Operator inference parameters
params.modelform = 'LQI';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

%% collect data for a series of trajectories with random inputs
num_inputs = 10;

if type == 1
    U_rand = rand(K,num_inputs);
elseif type == 2
    U_rand = 0.2*rand(K,num_inputs/2)-0.1;
    U_rand = [U_rand, zeros(K,num_inputs/2)];
else
    U_rand = rand(K,num_inputs)-0.5;
end

x_all = cell(num_inputs,1);
xdot_all = cell(num_inputs,1);
for i = 1:num_inputs
    if type == 1 || type == 2
        s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,i), IC);
    else
        IC_ = fic(ic_a(i));
        s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,i), IC_);
    end

    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:20;
err_inf = zeros(length(r_vals),1);  % relative state error for inferred model
err_int = zeros(length(r_vals),1);  % for intrusive model

% intrusive
rmax = max(r_vals);
Vr = U_svd(:,1:rmax);
Aint = Vr' * A * Vr;
Bint = Vr' * B;
Ln = elimat(N); Dr = dupmat(max(r_vals));
Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;

% op-inf (with stability check)
while true
    [operators] = inferOperators(X, U, Vr, params, R);
    Ahat = operators.A;
    Fhat = operators.F;
    Bhat = operators.B;

    % Check if the inferred operator is stable 
    lambda = eig(Ahat);
    Re_lambda = real(lambda);
    if all(Re_lambda(:) < 0)
        break;
    else
        warning("For mu = %f, order of r = %d is unstable. Decrementing max order.\n", mu, rmax);
        rmax = rmax - 1;
        Vr = U_svd(:,1:rmax);
    end
end

for j = 1:rmax
    r = r_vals(j);
    Vr = U_svd(:,1:r);

    Fhat_extract = extractF(Fhat, r);
    s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref,Vr'*IC);
    s_rec = Vr*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
    
    Fint_extract = extractF(Fint, r);
    s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref,Vr'*IC);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
end

%% Plotting
figure(3); clf
semilogy(r_vals(1:rmax),err_inf(1:rmax), DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals(1:rmax),err_int(1:rmax), DisplayName="int"); 
hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title("Burgers inferred model error, $\mu$ = "+num2str(mu),'Interpreter','LaTeX')
