function ECR = quadEnergyRate(A,X)
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