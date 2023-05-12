function xBu = controlEnergyRate(B,X,U)
   [~, K] = size(X);

    xBu = zeros(K,1);
    xBu(1) = X(:,1)' * B * 0;
    for i = 2:K
        xBu(i) = X(:,i)' * B * U(i-1);
    end
end