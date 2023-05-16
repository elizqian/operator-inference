function [CR,mmt] = constraintResidual_H(H)
    N = size(H,1);
    CR = 0;
    mmt = 0;
    for i = 1:N
        for j = 1:N
            for k = 1:N
                foo = H(i,N*(k-1)+j) + H(j,N*(k-1)+i) + H(k,N*(i-1)+j);
                CR = CR + abs(foo);
                mmt = mmt + foo;
            end
        end
    end
end