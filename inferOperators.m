function [operators] = inferOperators(X, U, Vr, params, dXdt)
% Infers linear, quadratic, bilinear, input, and constant matrix operators
% for state data in X projected onto basis in Vr and optional input data U.
% 
% INPUTS
% X         N-by-K full state data matrix
% U         K-by-m input data matrix
% Vr        N-by-r basis in which to learn reduced model
% params    struct with operator inference parameters - see PARAMS below
% rhs       N-b-K optional user-specified RHS for least-squares solve to be used
%           e.g., if user has Xdot data in the continuous time setting or
%           if data in X is non-uniformly spaced in time
%
% PARAMS
%   modelform   string indicating which terms to learn: e.g. 'LI' corresponds 
%               to linear model wit input; dxdt = A*x + B*u(t)
%               Options: 'L'inear, 'Q'uadratic, 'B'ilinear, 'I'nput, 'C'onstant;
%   modeltime   'continuous' e.g., if model form is dxdt = A*x OR 
%               'discrete' e.g., if model form is x(k+1) = A*x(k)
%   dt          timestep used to calculate state time deriv for
%               continuous-time models (use 0 if not needed)
%   lambda      L2 penalty weighting
%   scale       if true, scale data matrix to within [-1,1] before LS solve
%   implicit    true if data comes from implicit integrator
%
% OUTPUT
% operators     struct with inferred operators A, H, N, B, C. Terms that
%               are not part of the model are returned as empty matrices.

modelform = params.modelform;   
modeltime = params.modeltime;
dt        = params.dt;
lambda    = params.lambda;
scale     = params.scale;
implicit  = params.implicit;

m = size(U,2);

% if no right-hand side for LS problem is provided, calculate the rhs from
% state data based on specified model time
if nargin < 5
    K = size(X,2) - 1;
    switch modeltime
        case 'discrete'             
            rhs = X(:,2:end)'*Vr;
        case 'continuous'           
            Xdot = (X(:,2:end) - X(:,1:end-1))/dt;
            rhs = Xdot'*Vr;
    end
else
    K = size(X,2);
    rhs = dXdt'*Vr;
end

% get least-squares data matrix based on desired model form
[D,l,c,s,mr] = getDataMatrix(X,Vr,U,K,implicit,modelform);

% scale data before LS solve if desired
if scale
    scl = max(abs(D),[],1);
else
    scl = ones(1,size(D,2));
end
Dscl = D./scl;

% Solve LS problem and pull operators from result
temp = tikhonov(rhs,Dscl,lambda)';
temp = temp./scl;  % un-scale solution
operators.A = temp(:,1:l);
operators.F = temp(:,l+1:l+s);
operators.H = F2H(operators.F);
operators.N = temp(:,l+s+1:l+s+mr);
operators.B = temp(:,l+s+mr+1:l+s+mr+m);
operators.C = temp(:,l+s+mr+m+c);
end

%% builds data matrix based on desired form of learned model
function [D,l,c,s,mr] = getDataMatrix(X,Vr,U,K,implicit,modelform)
r = size(Vr,2);
if ~implicit
    Xhat = Vr'*X(:,1:K);
else
    Xhat = Vr'*X(:,2:K+1);
end

% if rhs contains B*u(t) input term
if contains(modelform,'I')
    U0 = U(1:K,:);
else
    U0 = [];
end

% if rhs contains quadratic H*kron(x,x) term
if contains(modelform,'Q')
    s = r*(r+1)/2;
    Xsq = get_x_sq(Xhat');
else
    s = 0;
    Xsq = [];
end

% if rhs contains constant term
if contains(modelform,'C')
    Y = ones(K,1);
    c = 1;
else
    Y = [];
    c = 0;
end

% if rhs contains linear A*x term
if ~contains(modelform,'L')
    Xhat = [];
    l = 0;
else
    l = r;
end

% if rhs contains bilinear N*x*u term
XU = [];
if contains(modelform,'B')
    for i = 1:m
        XU = [XU, Xhat'.*U0(:,i)];
    end
    mr = m*r;
else
    mr = 0;
end

D = [Xhat', Xsq, XU, U0, Y];
end