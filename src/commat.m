function K = commat(m, n)
% Produce commutation matrix from [Magnus and Neudecker, 1980]. The
% commutation matrix is defined to perform the relationship:
%               K * vec(A) = vec(A')  .
%
% INPUTS
% m         row dimension
% n         col dimension     
% 
% OUTPUTS
% K         commutation matrix with dim (m, n)
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023
%
% REFERENCE
% J. R. Magnus and H. Neudecker, “The Elimination Matrix: Some Lemmas and 
% Applications," SIAM. J. on Algebraic and Discrete Methods, vol. 1, no. 4,
% pp. 422–449, Dec. 1980, doi: 10.1137/0601049.

    mn = m*n;
    % determine permutation applied by K
    A = reshape(1:mn, m, n);
    v = reshape(A', 1, []);
    
    % apply this permutation to the rows (i.e. to each column) of identity matrix
    K = speye(mn);
    K = K(v,:);
end