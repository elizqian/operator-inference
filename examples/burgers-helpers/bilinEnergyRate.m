function xNxu = bilinEnergyRate(N,X,U)
   [~, K] = size(X);

    xNxu = zeros(K,1);
    xNxu(1) = X(:,1)' * N * X(:,1) * 0.0;
    for i = 2:K
        xNxu(i) = X(:,i)' * N * X(:,i) * U(i-1);
    end
end