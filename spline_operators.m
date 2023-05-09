function [A_interp, B_interp, F_interp] = spline_operators(ops, n, mus, mu_new)
    % Dimensions
    M = length(ops);
    l = size(ops{1}.B,2);
    s = n * (n+1) / 2;

    % Unpack matrices 
    A_all = zeros(n,n,M);
    B_all = zeros(n,l,M);
    F_all = zeros(n,s,M);
    for i = 1:M
        A_all(:,:,i) = ops{i}.A(1:n,1:n);
        B_all(:,:,i) = ops{i}.B(1:n,1:l);
        F_all(:,:,i) = ops{i}.F(1:n,1:s);
    end
    A_interp = spline(mus, A_all(:,:,:), mu_new);
    B_interp = spline(mus, B_all(:,:,:), mu_new);
    F_interp = spline(mus, F_all(:,:,:), mu_new);
end