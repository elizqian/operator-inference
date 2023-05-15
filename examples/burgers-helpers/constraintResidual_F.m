function CR = constraintResidual_F(F)
    N = size(F,1);
    CR = 0;
    for i = 1:N
        for j = 1:N
            for k = 1:N
                CR = (CR + abs(delta(j,k)*F(i,fidx(N,j,k)) ...
                    + delta(i,k)*F(j,fidx(N,i,k)) + delta(j,i)*F(k,fidx(N,j,i))));
            end
        end
    end
end

function del = delta(i,j)
    if i == j
        del = 1.0;
    else
        del = 0.5;
    end
end

function idx = fidx(n,j,k)
    if j >= k
        idx = (n - k/2)*(k - 1) + j;
    else
        idx = (n - j/2)*(j - 1) + k;
    end
end
