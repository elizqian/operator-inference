function K = commat(m, n)
    mn = m*n;
    % determine permutation applied by K
    A = reshape(1:mn, m, n);
    v = reshape(A', 1, []);
    
    % apply this permutation to the rows (i.e. to each column) of identity matrix
    K = speye(mn);
    K = K(v,:);
end