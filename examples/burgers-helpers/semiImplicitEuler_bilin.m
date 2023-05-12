%% semi-implicit Euler scheme for integrating learned model from zero initial condition
function s_hat = semiImplicitEuler_bilin(Ahat, Fhat, Bhat, Nhat, dt, u_input, IC)
    K = size(u_input,1);
    r = size(Ahat,1);
    s_hat = zeros(r,K+1); 
    ImdtA = eye(r) - dt*Ahat;

    s_hat(:,1) = IC;

    for i = 1:K
        ssq = get_x_sq(s_hat(:,i)')';
        s_hat(:,i+1) = ImdtA\(s_hat(:,i) + dt*Fhat*ssq + dt*Bhat*u_input(i) + dt*(Nhat*s_hat(:,i))*u_input(i));
        if any(isnan(s_hat(:,i)))
            warning(['ROM unstable at ',num2str(i-1),'th timestep'])
            break
        end
    end
end