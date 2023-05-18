%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler_noctrl(Ahat, Fhat, dt, K, IC)
    r = size(Ahat,1);
    s_hat = zeros(r,K); % initial state is zeros everywhere
    
    s_hat(:,1) = IC;

    ImdtA = eye(r) - dt*Ahat;
    for i = 2:K
        ssq = get_x_sq(s_hat(:,i-1)')';
        s_hat(:,i) = ImdtA\(s_hat(:,i-1) + dt*Fhat*ssq);
        if any(isnan(s_hat(:,i)))
            warning(['ROM unstable at ',num2str(i-1),'th timestep'])
            break
        end
    end
end