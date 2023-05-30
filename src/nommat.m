function N = nommat(m, n)
% Produce symmetric commutation matrix from [Magnus and Neudecker, 1980]. 
% This matrix is defined as:
%                   N = 0.5 * (I + K)  
% where I is the identity and K is the commutation matrix.
%
% INPUTS
% m         row dimension
% n         col dimension     
% 
% OUTPUTS
% N         symmetric commutation matrix with dim (m, n)
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023
%
% REFERENCE
% J. R. Magnus and H. Neudecker, “The Elimination Matrix: Some Lemmas and 
% Applications," SIAM. J. on Algebraic and Discrete Methods, vol. 1, no. 4,
% pp. 422–449, Dec. 1980, doi: 10.1137/0601049.

    mn = m*n;
    N = 0.5 * (speye(mn) + commat(m, n));
end