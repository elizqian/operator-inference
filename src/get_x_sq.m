function X2 = get_x_sq(X)
% takes data in X where each column is one variable and generates matrix X,
% where each column contains a quadratic term of the variables in X - this
% is different from a Kronecker product because there is no redundancy
%
% INPUT
% X     N-by-w data matrix where N is number of data points, w is number of
%       variables
%
% OUTPUT
% X2    N-by-(w*(w+1)/2) matrix of quadratic terms of data in X
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

w = size(X,2);
Xi = cell(w,1);
for j = 1:w
    Xi{j} = repmat(X(:,j),1,w-j+1).*X(:,j:w);
end
X2 = cat(2,Xi{:});
