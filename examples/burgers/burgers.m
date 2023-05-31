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
addpath('../../src/',"burgers-helpers/");

%% Problem set-up
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mus = 0.1:0.1:1.0; % diffusion coefficient
u_ref = ones(K,1);
IC = zeros(N,1);
Mp = 10;  % number of random inputs

% POD basis size and operators
r_vals = 1:15;
err_inf = zeros(length(mus),length(r_vals));
err_int = zeros(length(mus),length(r_vals));

% Store values
s_ref_all = cell(Mp,1);
infop_all = cell(Mp,1);
intop_all = cell(Mp,1);
Usvd_all = cell(10,1);

% Down-sampling
DS = 1;

%% Operator inference parameters
params.modelform = 'LQI';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

%% LEARN AND ANALYZE TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmin = max(r_vals);
for i = 1:length(mus)
    mu = mus(i);
    [A,B,F] = getBurgers_ABF_Matrices(N,1/(N-1),dt,mu);
    s_ref = semiImplicitEuler(A,F,B,dt,u_ref,IC);
    s_ref_all{i} = s_ref;

    %% collect data for a series of trajectories with random inputs
    U_rand = rand(K,Mp);
    x_all = cell(Mp,1);
    xdot_all = cell(Mp,1);
    for k = 1:Mp
        s_rand = semiImplicitEuler(A,F,B,dt,U_rand(:,k),IC);
        foo = s_rand(:,2:end);
        x_all{k} = foo(:,1:DS:end);  % down-sample and store
        bar = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
        xdot_all{k} = bar(:,1:DS:end);  % down-sample and store
    end
    
    X = cat(2,x_all{:});  % concatenate data from random trajectories
    R = cat(2,xdot_all{:}); 
    U = U_rand(1:DS:end,:);  % down-sample
    U = U(:);  % vectorize
    [U_svd,~,~] = svd(X,'econ');  % take SVD for POD basis

    % Intrusive operators
    rmax = max(r_vals);
    Vr = U_svd(:,1:rmax);
    Aint = Vr' * A * Vr;
    Bint = Vr' * B;
    Ln = elimat(N);  % elimination matrix of dim = N
    Dr = dupmat(max(r_vals));  % duplication matrix of dim = rmax
    Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;
    op_int.A = Aint; 
    op_int.B = Bint; 
    op_int.F = Fint;
    intop_all{i} = op_int;  % store the intrusive operator
    Usvd_all{i} = Vr;  % store the POD basis
    
    % Inferred operators with stability check
    while true
        [operators] = inferOperators(X, U, Vr, params, R);
        Ahat = operators.A;
        Fhat = operators.F;
        Bhat = operators.B;

        % Check if the inferred operator is stable 
        lambda = eig(Ahat);
        Re_lambda = real(lambda);
        if all(Re_lambda(:) < 0)
            infop_all{i} = operators;  % store operators
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
    
    % Update minimum stable reduced dim for plotting
    if rmin > rmax
        rmin = rmax;
    end
end

%% Plot relative state error
err_inf_avg = median(err_inf(:,1:rmin));
err_int_avg = median(err_int(:,1:rmin));

figure(1); clf;
semilogy(r_vals(1:rmin),err_inf_avg, DisplayName="opinf", Marker="o", MarkerSize=8); 
grid on; grid minor; hold on;
semilogy(r_vals(1:rmin),err_int_avg, DisplayName="int", Marker="x", MarkerSize=5); 
hold off; legend(Location="southwest");
xlabel('reduced model dimension $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error (Training)','Interpreter','LaTeX')

%% Verify Full Models
figure(2);
t = tiledlayout(2,Mp/2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for i = 1:Mp
    nexttile
    s = surf(linspace(0.0,T_end,K+1),linspace(0.0,1.0,N),s_ref_all{i}, ...
        'FaceAlpha',0.8);
    s.EdgeColor = 'none';
    xlabel("t");
    ylabel("x");
    title("$\mu$ = "+mus(i),Interpreter="latex")
end
