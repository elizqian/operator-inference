function N = nommat(m, n)
    mn = m*n;
    N = 0.5 * (speye(mn) + commat(m, n));
end