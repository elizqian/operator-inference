function Uinterp = spline_basis(U, mus, mu_new)
    % Dimensions
    M = length(U);
    [n,m] = size(U{1});

    % Unpack matrices 
    Uall = zeros(n,m,M); 
    for i = 1:M
        Uall(:,:,i) = U{i};
    end
    Uinterp = spline(mus, Uall(:,:,:), mu_new);
end