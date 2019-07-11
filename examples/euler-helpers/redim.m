% re-dimensionalizes non-dimensional Euler state data

% INPUTS
% state     non-dimensional state data
% ndc       struct for non-dim scaling constants {rho, u, l, t}
% type      'cons' or 'spv' indicates whether data in state is in 
%           conservative or specific volume representation

% OUTPUT
% state     dimensional state data

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 12 June 2019

function state = redim(state,ndc,type)
ntemp = size(state,1);
N = ntemp/3;

switch type
    case 'cons'
        state(1:N,:) = state(1:N,:)*ndc.rho;
        state(N+1:2*N,:) = state(N+1:2*N,:)*(ndc.rho*ndc.u);
        state(2*N+1:3*N,:) = state(2*N+1:3*N,:)*(ndc.rho*ndc.u^2);
    case 'spv'
        state(1:N,:) = state(1:N,:)*ndc.u;
        state(N+1:2*N,:) = state(N+1:2*N,:)*(ndc.rho*ndc.u^2);
        state(2*N+1:3*N,:) = state(2*N+1:3*N,:)/ndc.rho;
end