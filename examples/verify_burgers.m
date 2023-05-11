clear; clc;
addpath('../')

%% Problem set-up
N       = 2^7+1;        % num grid points
dt      = 1e-4;         % timestep
T_end   = 1;            % final time
K       = T_end/dt;     % num time steps

mu = 0.5;               % diffusion coefficient
type = 1;

% run FOM with input 1s to get reference trajectory
if type == 1
    u_ref = ones(K,1);
    IC = zeros(N,1);
else
    u_ref = zeros(K,1);
    IC = sin(pi * linspace(0,1,N))';
end
[A, B, F] = getBurgersMatrices(N,1/(N-1),dt,mu);
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

%% Slice of surface plot (time)
figure(2);
plot(0,0,MarkerSize=0.01,HandleVisibility="off")
hold on; grid on; grid minor; box on;
omega = linspace(0,1.0,N);
cmap = jet(length(1:floor(K/10):K+1));
ct = 1;
for i = 1:floor(K/10):K+1
    plot(omega,s_ref(:,i),Color=cmap(ct,:),DisplayName="$t="+num2str(i)+"$");
    ct = ct + 1;
end
hold off; legend(Interpreter="latex");
xlabel("\omega, space")
ylabel("x, velocity")
title("Burgers' plot sliced by time")

%% Slice of surface plot (space)
figure(3);
plot(0,0,MarkerSize=0.01,HandleVisibility="off")
hold on; grid on; grid minor; box on;
t = linspace(0,1.0,K+1);
cmap = jet(length(1:floor(N/10):N+1));
ct = 1;
for i = 1:floor(N/10):N+1
    plot(t,s_ref(i,:),Color=cmap(ct,:),DisplayName="$x="+num2str(i)+"$");
    ct = ct + 1;
end
hold off; legend(Interpreter="latex");
xlabel("t, time")
ylabel("x, velocity")
title("Burgers' plot sliced by space")

%% Plot the Energy
figure(2);
plot(linspace(0.0,T_end,K+1), vecnorm(s_ref))
xlabel("t, time")
ylabel("Energy")
grid on; grid minor; box on;
title("Energy over time of Burgers' Equation")

%% Check the Constraint Residual (H)
CR_h = crh(H);

%% Check the Constraint Residual (F)
CR_f = crf(F);

%% Plot the energy constraint residual (H)
ECR_h = ecr(H, s_ref);
ECR_f = ecr(F, s_ref);

figure(3);
plot(linspace(0.0,T_end,K+1), ECR_h, DisplayName="H", LineWidth=4, Color="b")
hold on; grid on; grid minor; box on;
plot(linspace(0.0,T_end,K+1), ECR_f, DisplayName="F", LineStyle="--", LineWidth=2, Color="g")
hold off; legend(Location="best");
xlabel("t, time")
ylabel("Energy")
title("Energy Constraint Residual over time of Burgers' Equation")

%% collect data for a series of trajectories with random inputs
num_inputs = 10;
U_rand = rand(K,num_inputs);
x_all = cell(num_inputs,1);
xdot_all = cell(num_inputs,1);
for i = 1:num_inputs
    s_rand = semiImplicitEuler(A, F, B, dt, U_rand(:,i), IC);
    x_all{i}    = s_rand(:,2:end);
    xdot_all{i} = (s_rand(:,2:end)-s_rand(:,1:end-1))/dt;
end

X = cat(2,x_all{:});        % concatenate data from random trajectories
R = cat(2,xdot_all{:});    
U = reshape(U_rand(:,1:num_inputs),K*num_inputs,1);

[U_svd,s_svd,~] = svd(X,'econ'); % take SVD for POD basis

%% for different basis sizes r, compute basis, learn model, and calculate state error 
r_vals = 1:25;
err_int = zeros(length(r_vals),1);  % for intrusive model
s_int_all = cell(length(r_vals));
s_tmp_all = cell(length(r_vals));

% intrusive
Vr = U_svd(:,1:max(r_vals));
Aint = Vr' * A * Vr;
Bint = Vr' * B;
Ln = elimat(N); Dr = dupmat(max(r_vals));
Fint = Vr' * F * Ln * kron(Vr,Vr) * Dr;

for j = 1:length(r_vals)
    r = r_vals(j);
    Vr = U_svd(:,1:r);
    
    Fint_extract = extractF(Fint, r);
    s_int = semiImplicitEuler(Aint(1:r,1:r),Fint_extract,Bint(1:r,:),dt,u_ref,Vr'*IC);
    s_tmp = Vr*s_int;
    err_int(j) = norm(s_tmp-s_ref,'fro')/norm(s_ref,'fro');
    s_int_all(i) = s_int;
    s_tmp_all(i) = s_tmp;
end

%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler(Ahat, Fhat, Bhat, dt, u_input, IC)
    K = size(u_input,1);
    r = size(Ahat,1);
    s_hat = zeros(r,K+1); % initial state is zeros everywhere
    
    s_hat(:,1) = IC;

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

%% Other Functions
function del = delta(i,j)
    if i == j
        del = 1.0;
    else
        del = 2.0;
    end
end

function idx = fidx(n,j,k)
    if j >= k
        idx = (n - k/2)*(k - 1) + j;
    else
        idx = (n - j/2)*(j - 1) + k;
    end
end

function CR = crf(F)
    N = size(F,1);
    CR = 0;
    for i = 1:N
        for j = 1:N
            for k = 1:N
                CR = (CR + delta(j,k)*F(i,fidx(N,j,k)) ...
                    + delta(i,k)*F(j,fidx(N,i,k)) + delta(j,i)*F(k,fidx(N,j,i)));
            end
        end
    end
end

function CR = crh(H)
    N = size(H,1);
    CR = 0;
    for i = 1:N
        for j = 1:N
            for k = 1:N
                CR = CR + H(i,N*(k-1)+j) + H(j,N*(k-1)+i) + H(k,N*(i-1)+j);
            end
        end
    end
end

function ECR = ecr(A,X)
    [n, K] = size(X);
    s = size(A,2);
    if s == n*(n+1)/2
        f = @(x) x' * A * get_x_sq(x')';
    elseif s == n^2
        f = @(x) x' * A * kron(x,x);
    else
        error("Unappropriate dimension for input matrix.");
    end

    ECR = zeros(K,1);
    for i = 1:K
        ECR(i) = f(X(:,i));
    end
end