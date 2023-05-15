% Elizabeth Qian (elizqian@mit.edu) 11 July 2019
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023 [EDITED]
% 
% THIS IS AN EXPERIMENTAL FILE FOR BILINEAR BURGERS EQUATION.
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
addpath('../',"burgers-helpers/");

%% Problem set-up
n       = 2^7+1;                  % num grid points
dt      = 1e-4;                   % timestep
T_end   = 1;                      % final time
K       = T_end/dt;               % num time steps
tspan = linspace(0.0,T_end,K+1);  % time span
sspan = linspace(0,1.0,n);        % spatial span

mu = 0.3;               % diffusion coefficient

type = 2;
% run FOM with input 1s to get reference trajectory
if type == 1
    u_ref = ones(K,1);
    IC = zeros(n,1);
else
    u_ref = zeros(K,1);
    IC = sin(pi * linspace(0,1,n))';
end

[A,B,F,N] = getBurgers_ABFN_Matrices(n,1/(n-1),mu);
H = F2Hs(F);
s_ref = semiImplicitEuler_bilin(A,F,B,N,dt,u_ref,IC);

%% Surface Plot for verification
figure(1);
s = surf(linspace(0.0,T_end,K+1),linspace(0.0,1.0,n),s_ref,'FaceAlpha',0.8,DisplayName="\mu="+num2str(mu));
s.EdgeColor = 'none';
xlabel("t, time");
ylabel("\omega, space");
zlabel("x(\omega,t), velocity")
axis tight
view(-73.25,38.649)
grid on
text(0.1,0.8,1,"\mu = "+num2str(mu),'FontSize',14);

%% Slice of surface plot (time)
figure(2);
plot(0,0,MarkerSize=0.01,HandleVisibility="off")
hold on; grid on; grid minor; box on;
omega = linspace(0,1.0,n);
cmap = jet(length(1:floor(K/10):K+1));
ct = 1;
for i = 1:floor(K/10):K+1
    plot(omega,s_ref(:,i),Color=cmap(ct,:),DisplayName="$t="+num2str(i)+"$");
    ct = ct + 1;
end
hold off; legend(Interpreter="latex");
xlabel("\omega, space")
ylabel("x, velocity")
title("Burgers' plot sliced by time \mu="+num2str(mu))

%% Slice of surface plot (space)
figure(3);
plot(0,0,MarkerSize=0.01,HandleVisibility="off")
hold on; grid on; grid minor; box on;
t = linspace(0,1.0,K+1);
cmap = jet(length(1:floor(n/10):n+1));
ct = 1;
for i = 1:floor(n/10):n+1
    plot(t,s_ref(i,:),Color=cmap(ct,:),DisplayName="$x="+num2str(i)+"$");
    ct = ct + 1;
end
hold off; legend(Interpreter="latex");
xlabel("t, time")
ylabel("x, velocity")
title("Burgers' plot sliced by space \mu="+num2str(mu))

%% Plot the Energy
figure(4);
plot(linspace(0.0,T_end,K+1), vecnorm(s_ref))
xlabel("t, time")
ylabel("Energy")
grid on; grid minor; box on;
title("Energy over time of Burgers' Equation")

%% Check the Constraint Residual (H)
CR_h = constraintResidual_H(H);

%% Check the Constraint Residual (F)
CR_f = constraintResidual_F(F);

%% Plot the Energy Rates
QER_h = quadEnergyRate(H, s_ref);
QER_f = quadEnergyRate(F, s_ref);
LER = linEnergyRate(A, s_ref);
CER = controlEnergyRate(B, s_ref, u_ref);
BER = bilinEnergyRate(N, s_ref, u_ref);

fig5 = figure(5);
fig5.Position = [250 500 2000 480];
t = tiledlayout(1,4,"TileSpacing","compact","Padding","compact");
nexttile;
    plot(tspan, QER_h, DisplayName="H", LineWidth=4, Color="b")
    hold on; grid on; grid minor; box on;
    plot(tspan, QER_f, DisplayName="F", LineStyle="--", LineWidth=2, Color="g")
    hold off; legend(Location="best");
    xlabel("t, time")
    ylabel("Energy Rate")
    title("Quadratic Energy Rate")
nexttile;
    plot(tspan, LER, Color="r", LineStyle=":", LineWidth=2, DisplayName="A")
    xlabel("t, time")
    ylabel("Energy Rate")
    title("Linear Energy Rate")
    grid on; grid minor; box on; legend(Location="best")
nexttile;
    plot(tspan, CER, Color="k", LineStyle="-.", LineWidth=2, DisplayName="B")
    xlabel("t, time")
    ylabel("Energy Rate")
    title("Control Energy Rate")
    grid on; grid minor; box on; legend(Location="best")
nexttile;
    plot(tspan, BER, Color="m", LineStyle="-.", LineWidth=2, DisplayName="N")
    xlabel("t, time")
    ylabel("Energy Rate")
    title("Bilinear Energy Rate")
    grid on; grid minor; box on; legend(Location="best")

figure(6);
TER = QER_h+LER+CER+BER;
semilogy(tspan, TER, Color="b", LineWidth=2)
xlabel("t, time")
ylabel("Energy Rate")
title("Total Energy Rate")
grid on; grid minor; box on;

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
    s_rand = semiImplicitEuler_bilin(A,F,B,N,dt,U_rand(:,i),IC);
    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:20;
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
    s_hat = semiImplicitEuler_bilin(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),Nhat(1:r,1:r),dt,u_ref,Vr'*IC);
    s_rec = Vr*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
    
    Fint_extract = extractF(Fint, r);
    s_int = semiImplicitEuler_bilin(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),Nint(1:r,1:r),dt,u_ref,Vr'*IC);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
end


%% Plotting
figure(7); clf;
semilogy(r_vals,err_inf, DisplayName="opinf"); grid on; grid minor; hold on;
semilogy(r_vals,err_int, DisplayName="int"); 
hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title('Burgers inferred model error','Interpreter','LaTeX')
