function CR = constraintResidual_H(H)
    N = size(H,1);
    CR = 0;
    for i = 1:N
        for j = 1:i
            for k = 1:j
                CR = CR + H(i,N*(k-1)+j) + H(j,N*(k-1)+i) + H(k,N*(i-1)+j);
            end
        end
    end
end