% generates a random smooth periodic initial condition for periodic-euler
% FOM via spline interpolation of Gaussian-sampled densities and velocities

% INPUT
% xgrid     mesh on which to generate initial condition

% OUTPUT
% init      struct that is an input to FOM.m

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function init = rand_init(xgrid)    % generate random smooth initial data
    N = length(xgrid);
    
    nip = 3; % number of interpolation points
    xx = downsample(xgrid,ceil(N/nip));
    
    % set initial rho to be cubic spline interp'd at random vals around 22
    yrho = 22 + randn(nip,1);   
    poly_rho = csape([xx; 2], [yrho; yrho(1)],'periodic');
    init.rho = ppval(poly_rho,xgrid);
    
    % set initial u as cubic spline interp'd at random vals around 100
    yu = 100+5*randn(nip,1);
    poly_u = csape([xx; 2],[yu; yu(1)],'periodic');
    init.u = ppval(poly_u,xgrid);
    
    % set initial p to be 1 bar everywhere
    init.p = 1e5*ones(N,1);
end