clear; clc; close all;
addpath('../', "burgers-helpers/");

%% Problem set-up
N       = 2^7+1;                  % num grid points
dt      = 1e-4;                   % timestep
T_end   = 1;                      % final time
K       = T_end/dt;               % num time steps
tspan = linspace(0.0,T_end,K+1);  % time span
sspan = linspace(0,1.0,N);        % spatial span

mu = 0.1;               % diffusion coefficient

type = 4;
% run FOM with input 1s to get reference trajectory
if type == 1
    u_ref = ones(K,1);
    IC = zeros(N,1);
elseif type == 2
    u_ref = zeros(K,1);
    IC = -exp(linspace(0,1,N) - 1)' .* sin(pi * linspace(0,1,N) - pi)';
elseif type == 3
    u_ref = zeros(K,1);
    IC = exp(-10*(3*linspace(0,1,N) - 1).^2)';
    fic = @(a) exp(-10*(a*linspace(0,1,N) - 1).^2);
    ic_a = linspace(2.2,4.0,10);

    [A, B, F] = getEPburgers_Matrices(N,1/(N-1),dt,mu);
    s_ref = semiImplicitEuler(A, F, B, dt, u_ref, IC);
elseif type == 4
%     IC = (sin(pi*linspace(0,1,N)/2 - pi/2).^2)';
%     IC = (pi/2*sin(2*pi*linspace(0,1,N)).^2 .* cos(pi*linspace(0,1,N)).^2)';
%     fic = @(a,b) sin(a*pi*linspace(0,1,N)/2 - pi/2*b).^2;
%     fic = @(a,b) a*pi/2*sin(b*2*pi*linspace(0,1,N)).^2 .* cos(pi*linspace(0,1,N)).^2;
%     ic_a = linspace(0.8,1.2,4);
%     ic_b = linspace(0.6,1.4,5);
%     [ic_a, ic_b] = meshgrid(linspace(0.8,1.2,5), linspace(0.7,1.3,5));
    
    IC = (sin(2*pi*linspace(0,1,N)))';
%     fic = @(a,b) a * sin(2*pi * linspace(0,1,N)) + b;
%     [ic_a, ic_b] = meshgrid([0.9 0.95 1.05 1.1], [-0.05 0.0 0.05]);

    ic_a = linspace(0.8,1.2,10);

    [A, F] = getEPburgers_Matrices(N,1/(N-1),mu);
    s_ref = semiImplicitEuler_noctrl(A, F, dt, K, IC);
end

% [A, B, F] = getBurgers_ABF_Matrices(N,1/(N-1),dt,mu);
% [A, B, F] = getEPburgers_Matrices(N,1/(N-1),dt,mu);
% [A, B, F] = getCFburgers_Matrices(N,1/(N-1),dt,mu);

H = F2Hs(F);

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

%% Slice of surface plot (time)
figure(2);
plot(0,0,MarkerSize=0.01,HandleVisibility="off")
hold on; grid on; grid minor; box on;
cmap = jet(length(1:floor(K/10):K+1));
ct = 1;
for i = 1:floor(K/10):K+1
    plot(sspan,s_ref(:,i),Color=cmap(ct,:),DisplayName="$t="+num2str(i)+"$");
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
cmap = jet(length(1:floor(N/10):N+1));
ct = 1;
for i = 1:floor(N/10):N+1
    plot(tspan,s_ref(i,:),Color=cmap(ct,:),DisplayName="$x="+num2str(i)+"$");
    ct = ct + 1;
end
hold off; legend(Interpreter="latex");
xlabel("t, time")
ylabel("x, velocity")
title("Burgers' plot sliced by space \mu="+num2str(mu))

%% Plot the Energy
figure(4);
plot(tspan, vecnorm(s_ref).^2/2)
xlabel("t, time")
ylabel("Energy")
grid on; grid minor; box on;
title("Energy over time of Burgers' Equation")

%% Check the Constraint Residual (H)
[CR.fom.H, mmt.fom.H] = constraintResidual_H(H);

%% Check the Constraint Residual (F)
[CR.fom.F, mmt.fom.F] = constraintResidual_F(F);

%% Elizabeth's Sandbox
get_h = @(i,j,k) H(i,(k-1)*N+j);
derp = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for k = 1:N
            derp(i,j,k) = get_h(i,j,k) + get_h(j,i,k) + get_h(k,j,i);
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

%% Plot the Energy Rates
QER_h = quadEnergyRate(H, s_ref);
QER_f = quadEnergyRate(F, s_ref);
LER = linEnergyRate(A, s_ref);
% CER = controlEnergyRate(B, s_ref, u_ref);

fig5 = figure(5);
fig5.Position = [500 500 1500 480];
t = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
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

% figure(6);
% TER = QER_h+LER+CER;
% semilogy(tspan, TER, Color="b", LineWidth=2)
% xlabel("t, time")
% ylabel("Energy Rate")
% title("Total Energy Rate")
% grid on; grid minor; box on;

%%
figure(199);
plot(tspan, QER_h, DisplayName="H", LineWidth=4, Color="b")
hold on; grid on; grid minor; box on;
plot(tspan, QER_f, DisplayName="F", LineStyle="--", LineWidth=2, Color="g")
hold off; legend(Location="best");
ylim([-8e-15, 6e-15])
xlabel("t, time")
ylabel("Energy Rate")
title("Quadratic Energy Rate")


%% collect data for a series of trajectories with random inputs
num_inputs = 10;
if type == 1
    U_rand = rand(K,num_inputs);
elseif type == 2
    U_rand = 0.2*rand(K,num_inputs/2)-0.1;
    U_rand = [U_rand, zeros(K,num_inputs/2)];
elseif type == 3
    U_rand = rand(K,num_inputs)-0.5;
elseif type == 4
    % no control
end

x_all = cell(num_inputs,1);
xdot_all = cell(num_inputs,1);
for i = 1:num_inputs
    if type == 1 || type == 2
        s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,i), IC);
    elseif type == 3
        IC_ = fic(ic_a(i));
        s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,i), IC_);
    elseif type == 4
        s_rand = semiImplicitEuler_noctrl(A, F, dt, K, ic_a(i) * IC);
    end

    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
% U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% Operator inference parameters
params.modelform = 'LQ';           % model is linear-quadratic with input term
params.modeltime = 'continuous';    % learn time-continuous model
params.dt        = dt;              % timestep to compute state time deriv
params.ddt_order = '1ex';           % explicit 1st order timestep scheme

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:20;
err_inf = zeros(length(r_vals),1);  % relative state error for inferred model
err_int = zeros(length(r_vals),1);  % for intrusive model

% intrusive
rmax = max(r_vals);
Vr = U_svd(:,1:rmax);
Aint = Vr' * A * Vr;
% Bint = Vr' * B;
Ln = elimat(N); Dr = dupmat(rmax);
Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;
Hint = Vr' * H * kron(Vr,Vr);

[operators] = inferOperators(X, NaN, Vr, params, R);
Ahat = operators.A;
Fhat = operators.F;

% op-inf (with stability check)
% while true
%     [operators] = inferOperators(X, NaN, Vr, params, R);
%     Ahat = operators.A;
%     Fhat = operators.F;
% %     Bhat = operators.B;
%     
%     % Check if the inferred operator is stable 
%     lambda = eig(Ahat);
%     Re_lambda = real(lambda);
%     if all(Re_lambda(:) < 0)
%         break;
%     else
%         warning("For mu = %f, order of r = %d is unstable. Decrementing max order.\n", mu, rmax);
%         rmax = rmax - 1;
%         Vr = U_svd(:,1:rmax);
%     end
% end

for j = 1:rmax
    r = r_vals(j);
    Vr = U_svd(:,1:r);

    Fhat_extract = extractF(Fhat, r);
%     s_hat = semiImplicitEuler(Ahat(1:r,1:r),Fhat_extract,Bhat(1:r,:),dt,u_ref,Vr'*IC);
    s_hat = semiImplicitEuler_noctrl(Ahat(1:r,1:r),Fhat_extract,dt,K,Vr'*IC);
    s_rec = Vr*s_hat;
    err_inf(j) = norm(s_rec-s_ref,'fro')/norm(s_ref,'fro');
    
    Fint_extract = extractF(Fint, r);
%     s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref,Vr'*IC);
    s_int = semiImplicitEuler_noctrl(Aint(1:r,1:r),Fint_extract,dt,K,Vr'*IC);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
end

%% Plotting
figure(7); clf
semilogy(r_vals(1:rmax),err_inf(1:rmax), DisplayName="opinf", Marker="o", MarkerSize=8); 
grid on; grid minor; hold on;
semilogy(r_vals(1:rmax),err_int(1:rmax), DisplayName="int", Marker="x", MarkerSize=5); 
hold off; legend(Location="southwest");
xlabel('Model size $r$','Interpreter','LaTeX')
ylabel('Relative state reconstruction error','Interpreter','LaTeX')
title("Burgers inferred model error, $\mu$ = "+num2str(mu),'Interpreter','LaTeX')

%% Check the Constraint Residual for inferred (H)
[CR.inf.H, mmt.inf.H] = constraintResidual_H(operators.H);

%% Check the Constraint Residual for inferred (F)
[CR.inf.F, mmt.inf.F] = constraintResidual_F(operators.F);

%% Plot the Energy Rates
QER_h = quadEnergyRate(operators.H, U_svd(:,1:rmax)' * s_ref);
QER_f = quadEnergyRate(operators.F, U_svd(:,1:rmax)' * s_ref);
LER = linEnergyRate(operators.A, U_svd(:,1:rmax)' * s_ref);
% CER = controlEnergyRate(operators.B, U_svd(:,1:rmax)' * s_ref, u_ref);

fig8 = figure(8);
fig8.Position = [500 500 1500 480];
t = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
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
%     plot(tspan, CER, Color="k", LineStyle="-.", LineWidth=2, DisplayName="B")
%     xlabel("t, time")
%     ylabel("Energy Rate")
%     title("Control Energy Rate")
%     grid on; grid minor; box on; legend(Location="best")

% figure(9);
% TER = QER_h+LER+CER;
% semilogy(tspan, TER, Color="b", LineWidth=2)
% xlabel("t, time")
% ylabel("Energy Rate")
% title("Total Energy Rate")
% grid on; grid minor; box on;

%%
figure(220);
semilogy(tspan, QER_h, DisplayName="H", LineWidth=4, Color="b")
hold on; grid on; grid minor; box on;
plot(tspan, QER_f, DisplayName="F", LineStyle="--", LineWidth=2, Color="g")
hold off; legend(Location="best");
xlabel("t, time")
ylabel("Energy Rate")
title("Quadratic Energy Rate")


%% Check the Constraint Residual for intrusive (H)
[CR.int.H, mmt.int.H] = constraintResidual_H(Hint);

%% Check the Constraint Residual for intrusive (F)
[CR.int.F, mmt.int.F] = constraintResidual_F(Fint);

%% Plot the Energy Rates
rmax = max(r_vals);
QER_h = quadEnergyRate(Hint, U_svd(:,1:rmax)' * s_ref);
QER_f = quadEnergyRate(Fint, U_svd(:,1:rmax)' * s_ref);
LER = linEnergyRate(Aint, U_svd(:,1:rmax)' * s_ref);
% CER = controlEnergyRate(Bint, U_svd(:,1:rmax)' * s_ref, u_ref);

fig10 = figure(10);
fig10.Position = [500 500 1500 480];
t = tiledlayout(1,3,"TileSpacing","compact","Padding","compact");
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

figure(11);
TER = QER_h+LER+CER;
semilogy(tspan, TER, Color="b", LineWidth=2)
xlabel("t, time")
ylabel("Energy Rate")
title("Total Energy Rate")
grid on; grid minor; box on;

