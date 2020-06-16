function [V, D, output] = gevd_bal(A, B, varargin)
%GEVD_BAL Generalized eigenvalue decomposition with balancing.
%   E = GEVD_BAL(A, B) computes the generalized eigenvalues of square matrices
%   A and B by first balancing the problem using the Skinkhorn-Knopp
%   algorithm as described in [1], i.e. it is attempted to scale
%   |A|.^2+|B|.^2 into a doubly stochastic matrix. 
%
%   [V, D] = GEVD_BAL(A, B) returns the generalized eigenvectors in V and the
%   generalized eigenvalues on the diagonal of D.
%
%   [V, D, output] = GEVD_BAL(A, B) also returns a struct with the termination
%   status of the optimization loop and contains the following fields:
%      iterations               - The number of iterations performed.
%      info                     - The reason for termination:
%                                 1 = relative improvement < TolBal
%                                 3 = maximum number of iterations reached
%
%   GEVD_BAL(A, B, options) where options is a struct, or contains is a series
%   of key-value pairs, can be used to set the following options:
%
%      MaxIter = 4*size(A, 1)   - The maximum number of iterations for
%                                 balancing. 
%      TolBal = 1e-6            - The balancing happens iteratively until the
%                                 relative improvement falls below TolBal.
% 
%   See also eig.
    
% References:
% [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Numerical solution of
%     bivariate and polyanalytic polynomial systems", SIAM Journal on
%     Numerical Analysis, Vol. 52, No. 3, 2014, pp. 1551-1572 .
%
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2014/12/02   NV      Initial version
    
    p = inputParser;
    p.addOptional('TolBal', 1e-6);
    p.addOptional('MaxIter', 4*size(A, 1));
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    if any(size(A) ~= size(B))
        error('gevd_bal:size', 'size(A) should be equal to size(B)');
    end
    if size(A,1) ~= size(A, 2)
        error('gevd_bal:size', 'A and B should be square');
    end
    
    M = A.*conj(A) + B.*conj(B);
    n = size(A, 1);
    dl = ones(1, n);
    dr = ones(n, 1);
    
    output.info = 3;
    
    t = (M*dr).';
    dy = max(abs(t-1)); % dl = 1
    for k = 1:options.MaxIter
        % compute updates
        dl = 1./t;
        t = dl*M;
        dr = 1./t.';
        t = (M*dr).';
        % check stopping criterion
        dy1 = dy;
        dy = max(abs(dl.*t-1));
        if 1-dy/dy1 <= options.TolBal, 
            output.info = 1;
            break; 
        end
    end
    output.iterations = k;
    
    % Balance problem 
    if ~any(isnan(dl)) && ~any(isnan(dr))
        dl = 2.^round(0.5*log2(abs(dl)));
        dr = 2.^round(0.5*log2(abs(dr)));
        A = bsxfun(@times,bsxfun(@times,A,dl.'),dr.');
        B = bsxfun(@times,bsxfun(@times,B,dl.'),dr.');
    end
    
    % compute generalized eigenvalues
    if nargout <= 1; 
        V = eig(A,B); 
    else
        [V, D] = eig(A, B);
        if ~any(isnan(dr))
            V = bsxfun(@times, V, dr);
        end
    end
end
