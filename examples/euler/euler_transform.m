% Elizabeth Qian (elizqian@mit.edu) 11 July 2019
% -----------------------------------------------------------------------
% Uses Transform & Learn approach to learn a quadratic reduced model for 
% the 1D Euler equations on a periodic domain. Data from a simulation in 
% the conservative variables are first transformed to the quadratic 
% specific volume representation and a quadratic model is learned. 
%
% Generation of random smooth initial conditions requires MATLAB Curve
% Fitting Toolbox
%
% Based on:
%   Qian, E., Kramer, B., Marques, A. and Willcox, K., "Transform & Learn:
%   A data-driven approach to nonlinear model reduction." In AIAA Aviation 
%   2019 Forum, June 17-21, Dallas, TX.
%
% See also:
%   Peherstorfer, B. and Willcox, K., "Data-driven operator inference for
%   non-intrusive projection-based model reduction." Computer Methods in
%   Applied Mechanics and Engineering, 306:196-215, 2016.

%% SETUP
clear
addpath('euler-helpers','../../src/')

% FOM parameters
N       = 200;
dt      = 1e-5;
Tfinal  = 1e-2;
xgrid   = linspace(0,2,N+1)';
xgrid   = xgrid(1:end-1);

% non-dimensionalization scaling constants
ndc.rho = 10;   
ndc.u = 100;
ndc.l = 2;
ndc.t = ndc.l/ndc.u;

%% ACQUIRE TRAINING DATA FROM NONLINEAR FOM
% run FOM to acquire training state data
M = 20;
s_train = cell(M,1);
for i = 1:M
    init = rand_init(xgrid);
    [s,init,time,runtime] = FOM(xgrid,dt,init,Tfinal);
    s_train{i} = [init.s0, s];
end

%% TRANSFORM DATA TO QUADRATIC STATE
% process data to get non-dim transformed state & time-derivative data
w_all = cell(M,1);
r_all = cell(M,1);
for i = 1:M
    % transform state to specific volume representation
    w        = nonlin2quad(s_train{i});     
    
    % non-dimensionalize state so that POD basis is better
    w_nd     = nondim(w,ndc,'spv');
    
    % compute non-dim residual in sp. vol. variables
    r_nd     = (w_nd(:,2:end)-w_nd(:,1:end-1))/(dt/ndc.t);
    
    w_all{i} = w_nd(:,1:end-1);
    r_all{i} = r_nd;
end
W = cat(2,w_all{:});
R = cat(2,r_all{:});

%% FORMULATE AND SOLVE LEARNING PROBLEM
% take SVD to get POD modes
[Usvd,s,~] = svd(W,'econ');

% set POD basis size (same as dimension of model to be learned)
r = 20;
Vr = Usvd(:,1:r);

params.modelform = 'Q';
params.lambda    = 0.001;

operators = inferOperators(W,[],Vr,params,R);
H         = operators.H;

disp(['Transform & Learn model of size r = ',num2str(r),' learned from ',num2str(M),' random training trajectories.'])
%% COMPUTE TRANFORM & LEARN TRAINING ERROR
ndtime = time/ndc.t;
f_rom = @(w) H*kron(w,w);     % define nondim ROM residual

training_err = zeros(M,1);
for i = 1:M
    what0 = Vr'*w_all{i}(:,1);                  % ROM initial condition
    w_rom = forwardEuler(f_rom,what0,ndtime);   % non-dimensional integration
    w_rec = Vr*w_rom;                           % reconstruct high-dim state
    
    % compute reference solution in nondim sp. vol. variables
    w_ref = nondim(nonlin2quad(s_train{i}(:,2:end)),ndc,'spv');
    
    training_err(i) = norm(w_rec - w_ref,'fro')/norm(w_ref,'fro');
end

disp(['Mean relative training error: ',num2str(mean(training_err))])

%% COMPUTE TEST DATA SET
s_test = cell(M,1);
for i = 1:M
    init = rand_init(xgrid);
    [s,init,time,runtime] = FOM(xgrid,dt,init,Tfinal);
    s_test{i} = [init.s0, s];
end

%% COMPUTE TEST ERROR
test_err = zeros(M,1);
for i = 1:M
    % ROM initial condition 
    what0 = Vr'*nondim(nonlin2quad(s_test{i}(:,1)),ndc,'spv');  
     
    w_rom = forwardEuler(f_rom,what0,ndtime);   % non-dimensional integration
    w_rec = Vr*w_rom;                           % reconstruct high-dim state
    
    w_ref = nondim(nonlin2quad(s_test{i}(:,2:end)),ndc,'spv');
    test_err(i) = norm(w_rec - w_ref,'fro')/norm(w_ref,'fro');
end

disp(['Mean relative test error:     ',num2str(mean(test_err))])

%% Forward Euler integrator
function y_all = forwardEuler(f,y0,t)
K = length(t)-1;
n = length(y0);

y_all = zeros(n,K);
y_all(:,1) = y0 + (t(2)-t(1))*f(y0);
for i = 2:K
    y_all(:,i) = y_all(:,i-1) + (t(i+1)-t(i))*f(y_all(:,i-1));
end
end
