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