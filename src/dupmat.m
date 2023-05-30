function D2 = dupmat(n)
% Produce duplication matrix from [Magnus and Neudecker, 1980]. The
% duplication matrix is defined to perform the relationship:
%               D * vech(A) = vec(A)  if  A = A'  .
%
% INPUTS
% n         row, col dimension     
% 
% OUTPUTS
% D         duplication matrix with dim (n, n)
%
% AUTHOR
% Tomoki Koike (tkoike3@gatech.edu) 12 May 2023 [Author of function file]
% Jan: https://www.mathworks.com/matlabcentral/answers/473737-efficient-
% algorithm-for-a-duplication-matrix#answer_385153 [Credit for algorithm]
%
% REFERENCE
% J. R. Magnus and H. Neudecker, “The Elimination Matrix: Some Lemmas and 
% Applications," SIAM. J. on Algebraic and Discrete Methods, vol. 1, no. 4,
% pp. 422–449, Dec. 1980, doi: 10.1137/0601049.

  m   = n * (n + 1) / 2;
  nsq = n^2;
  r   = 1;
  a   = 1;
  v   = zeros(1, nsq);
  cn  = cumsum(n:-1:2);  
  for i = 1:n
     v(r:r + i - 2) = i - n + cn(1:i - 1);   
     r = r + i - 1;
     
     v(r:r + n - i) = a:a + n - i;
     r = r + n - i + 1;
     a = a + n - i + 1;
  end
  
  D2 = sparse(1:nsq, v, 1, nsq, m);
end