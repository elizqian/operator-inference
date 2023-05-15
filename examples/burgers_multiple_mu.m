% Elizabeth Qian (elizqian@mit.edu) 11 July 2019
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023 [EDITED]
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


% Operator inference parameters
params.modelform = 'LQI';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

type = 2;
% FOM with 2 types of inputs and BCs
if type == 1
    mus = 0.1:0.1:1.0; % diffusion coefficient
    u_ref = ones(K,1);
    IC = zeros(N,1);
else
    % mus = 0.05:0.05:0.5; % diffusion coefficient
    mus = 0.1:0.1:1.0; % diffusion coefficient
    u_ref = zeros(K,1);
    IC = exp(-10*(5*linspace(0,1,N) - 1).^2)';
    fic = @(a) exp(-10*(a*linspace(0,1,N) - 1).^2);
    ic_a = linspace(3,7.5,10);
end

% Settings
M = length(mus);
r_vals = 1:20;
err_inf = zeros(length(mus),length(r_vals));
err_int = zeros(length(mus),length(r_vals));

% Store all s_ref values for verification
s_ref_all = cell(M,1);

% Store all operators to use for test data
infop_all = cell(M,1);
intop_all = cell(M,1);

%% LEARN AND ANALYZE TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin = max(r_vals);
for i = 1:M
    mu = mus(i);
    [A,B,F] = getBurgers_ABF_Matrices(N,1/(N-1),dt,mu);
    s_ref = semiImplicitEuler(A,F,B,dt,u_ref,IC);
    s_ref_all{i} = s_ref;

    %% collect data for a series of trajectories with random inputs
    num_inputs = 10;
    if type == 1
        U_rand = rand(K,num_inputs);
    elseif type == 2
        U_rand = rand(K,num_inputs)-0.5;
    end

    x_all = cell(num_inputs,1);
    xdot_all = cell(num_inputs,1);
    for k = 1:num_inputs
        if type == 1
            s_rand = semiImplicitEuler(A,F,B,dt,U_rand(:,k),IC);
        else
            IC_ = fic(ic_a(k));
            s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,k), IC_); 
        end
        x_all{k}    = s_rand(:,2:end);
        xdot_all{k} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
    end
    
    X = cat(2,x_all{:});  % concatenate data from random trajectories
    R = cat(2,xdot_all{:});    
    U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);
    [U_svd,~,~] = svd(X,'econ');  % take SVD for POD basis

    % intrusive
    rmax = max(r_vals);
    Vr = U_svd(:,1:rmax);
    Aint = Vr' * A * Vr;
    Bint = Vr' * B;
    Ln = elimat(N); Dr = dupmat(max(r_vals));
    Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;
    op_int.A = Aint; op_int.B = Bint; op_int.F = Fint;
    intop_all{i} = op_int;
    
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
            infop_all{i} = operators;
            break;
        else
            warning("For mu = %f, order of r = %d is unstable. Decrementing max order.\n", mu, rmax);
            rmax = rmax - 1;
            Vr = U_svd(:,1:rmax);
        end
    end

    %% For different basis sizes r, compute basis, learn model, and calculate state error 
    for j = 1:rmax
        r = r_vals(j);
        Vr = U_svd(:,1:r);
        
        % Extract operators for inferred model
        Fhat_extract = extractF(Fhat, r);
        s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref,Vr'*IC);
        s_rec = Vr*s_hat;
        err_inf(i,j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
        
        % Extract operators for intrusive model
        Fint_extract = extractF(Fint, r);
        s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref,Vr'*IC);
        s_tmp = Vr*s_int;
        err_int(i,j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
    end

    if rmin > rmax
        rmin = rmax;
    end
end

%% Plot relative state error
err_inf_avg = mean(err_inf(:,1:rmin),"omitnan");
err_int_avg = mean(err_int(:,1:rmin),"omitnan");

figure(1); clf;
semilogy(r_vals(1:rmin),err_inf_avg, DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals(1:rmin),err_int_avg, DisplayName="int"); hold off; legend(Location="southwest");
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


% %% ANALYZE TEST DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mus_new = rand(5,1);
% M_new = length(mus_new);
% 
% err_inf_new = zeros(M_new,length(r_vals));
% err_int_new = zeros(M_new,length(r_vals));
% 
% %%
% for i = 1:M_new
%     mu = mus_new(i);
%     [A,B,F] = getBurgersMatrices(N,1/(N-1),dt,mu);
%     s_ref = semiImplicitEuler(A,F,B,dt,u_ref);
% 
%     % Interpolate the intrusive and inferred operator for random mu value
%     [Aint, Bint, Fint] = spline_operators(intop_all, max(r_vals), mus, mu);
%     [Ainf, Binf, Finf] = spline_operators(infop_all, max(r_vals), mus, mu);
% 
%     %% For different basis sizes r, compute basis, learn model, and calculate state error 
%     for j = 1:length(r_vals)
%         r = r_vals(j);
%         Vr = U_svd(:,1:r);
%         
%         % Extract operators for inferred model
%         Fhat_extract = extractF(Fhat, r);
%         s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref);
%         s_rec = Vr*s_hat;
%         err_inf_new(i,j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
%         
%         % Extract operators for intrusive model
%         Fint_extract = extractF(Fint, r);
%         s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref);
%         s_tmp = Vr*s_int;
%         err_int_new(i,j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
%     end
% end
% 
% %% Plotting
% err_inf_avg_new = mean(err_inf_new);
% err_int_avg_new = mean(err_int_new);
% 
% figure(3); clf;
% semilogy(r_vals,err_inf_avg_new, DisplayName="opinf"); grid on; grid minor; hold on;
% semilogy(r_vals,err_int_avg_new, DisplayName="int"); hold off; legend(Location="southwest");
% xlabel('Model size $r$','Interpreter','LaTeX')
% ylabel('Relative state reconstruction error','Interpreter','LaTeX')
% title('Burgers inferred model error (Test)','Interpreter','LaTeX')
