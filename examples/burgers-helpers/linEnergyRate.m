function xAx = linEnergyRate(A,X)
   [~, K] = size(X);

    xAx = zeros(K,1);
    for i = 1:K
        xAx(i) = X(:,i)' * A * X(:,i);
    end
end