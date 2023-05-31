function H = F2Hs(F)
% converts F matrix operator which operates on get_x_sq output (no
% Kronecker redundancy) to tensorized matrix operator H which operates on
% kron(x,x) output (with redundancy) [Same as F2H.m but faster]
%
% INPUT
% F     w-by-(w*(w+1)/2) matrix operator that acts on get_x_sq.m output
%
% OUTPUT
% H     w-by-w^2 tensorized matrix operator that acts on kron(x,x)
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023
%
% REFERENCE
% J. R. Magnus and H. Neudecker, “The Elimination Matrix: Some Lemmas and 
% Applications," SIAM. J. on Algebraic and Discrete Methods, vol. 1, no. 4,
% pp. 422–449, Dec. 1980, doi: 10.1137/0601049.

    n = size(F,1);
    H = F * elimat(n) * nommat(n,n);
end