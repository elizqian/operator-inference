function L = elimat(m)
% Produce elimination matrix from [Magnus and Neudecker, 1980]. The
% elimination matrix is defined to perform the relationship:
%               L * vec(A) = vech(A)  if  A = A' .
%
% INPUTS
% m         row, col dimension
% 
% OUTPUTS
% L         elimination matrix with dim (m, m)
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023 [Author of function file]
% James Tursa: https://www.mathworks.com/matlabcentral/answers/274256-
% algorithm-for-an-elimination-matrix#answer_214142 [Credit for algorithm]
%
% REFERENCE
% J. R. Magnus and H. Neudecker, “The Elimination Matrix: Some Lemmas and 
% Applications," SIAM. J. on Algebraic and Discrete Methods, vol. 1, no. 4,
% pp. 422–449, Dec. 1980, doi: 10.1137/0601049.

  T = tril(ones(m)); % Lower triangle of 1's
  f = find(T(:)); % Get linear indexes of 1's
  k = m*(m+1)/2; % Row size of L
  m2 = m*m; % Colunm size of L
  L = zeros(m2,k); % Start with L'
  x = f + m2*(0:k-1)'; % Linear indexes of the 1's within L'
  L(x) = 1; % Put the 1's in place
  L = sparse(L'); % Now transpose to actual L
end