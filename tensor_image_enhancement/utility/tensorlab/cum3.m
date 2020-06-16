function c3 = cum3(X,prewhiten)
%CUM3 Third-order cumulant tensor.
%   c3 = cum3(X) computes the third-order cumulant of a matrix X in which
%   each row is an observation and each column is a variable. Herein,
%
%      c3(i,j,k) = E[xi.*conj(xj).*conj(xk)]
%
%   where the expectation E is approximated by the arithmetic mean and xi
%   is the i-th mean centered variable, X(:,i)-mean(X(:,i)) (and
%   analogously for xj, xk and xl).
%
%   cum3(X,'prewhiten') applies a linear transformation to the columns of X
%   so that the covariance matrix of the new matrix is the identity matrix
%   before computing its fourth-order cumulant.
%
%   See also cov, scov, cum4.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] P. McCullagh, "Tensor Methods in Statistics," Chapman and Hall,
%       London, 1987.
%   [2] C. Nikias, A. Petropulu, "Higher-Order Spectra Analysis: A 
%       Nonlinear Signal Processing Framework," Prentice Hall, 1993.

% Check the prewhiten option.
if nargin < 2, prewhiten = false; end
if ischar(prewhiten), prewhiten = strcmpi(prewhiten,'prewhiten'); end

% Center the variables.
X = bsxfun(@minus,X,mean(X,1));

% Apply a prewhitening to X if requested.
n = size(X,1);
if prewhiten
    [U,S,~] = svd(X,'econ');
    X = U*(S*pinv(S))*sqrt(n-1);
end

% Compute c3 = E[xi*conj(xj)*conj(xk)].
c3 = mean(bsxfun(@times,bsxfun(@times,permute(X,[2 3 4 1]), ...
    permute(conj(X),[3 2 4 1])),permute(conj(X),[3 4 2 1])),4);
