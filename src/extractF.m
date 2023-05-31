function Fhat = extractF(F,r)
% Extract the lower dimension (r) of the quadratic operator, F.
%
% INPUTS
% F         w-by-(w*(w+1)/2) matrix operator that acts on get_x_sq.m output
% r         reduced dimension
% 
% OUTPUTS
% Fhat      extracted quadratic operator
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023 

    N = size(F,1);
    xsq_idx = zeros(1,N);
    for n = 1:N
        xsq_idx(n) = 1+(N+1)*(n-1)-n*(n-1)/2;
    end
    extract_idx = [];
    for i = 1:r
        x = xsq_idx(i);
        extract_idx = [extract_idx, x:(x+r-i)];
    end
    Fhat = F(1:r, extract_idx);
end