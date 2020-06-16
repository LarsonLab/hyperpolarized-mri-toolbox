function T = cpdgen(U,varargin)
%CPDGEN Generate full tensor given a polyadic decomposition.
%   T = CPDGEN(U) computes the tensor T as the sum of R rank-one tensors
%   defined by the columns of the factor matrices U{n}.
%
%   T = CPDGEN(U,ind) computes only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = CPDGEN(U,i,j,k,...) computes only the elements T(i,j,k,...) from the
%   tensor. The number of indices should match the order of the tensor, given
%   by length(U). The colon operator can be used as a string, e.g., for a
%   third order tensor, cpdgen(U,i,':',k) can be used. 
%
%   See also btdgen, lmlragen, ll1gen, ttgen, cpdres.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/01/13  NV   Added incomplete cpdgen
%   - 2014/12/18  NV   Reordering of Khatri-Rao products to reduce memory
%                      footprint

size_tens = cellfun('size', U(:).', 1);
N = length(size_tens);
if nargin == 1
    % Split modes and balance them such that both kr-products have a similar size
    modes = cumprod(size_tens) <= round(sqrt(prod(size_tens)));
    modes(1) = 1;
    modes(end) = 0;
    % Compute tensor
    tmp = U(modes);
    l = kr(tmp(end:-1:1));
    tmp = U(~modes);
    r = kr(tmp(end:-1:1));
    res = l*r.';
    T = reshape(res, size_tens); 
elseif nargin == 2
    if ischar(varargin{1}) && strcmp(varargin{1},':')
        T = reshape(cpdgen(U), [], 1);
        return;
    end
    size_output = size(varargin{1});
    if nargin == 2 % linear indexing
        sub = cell(1, N);
        [sub{:}] = ind2sub(size_tens, varargin{1});
    end
    T = U{1}(sub{1},:);
    for n = 2:N
        T = T .* U{n}(sub{n},:);
    end
    T = T*ones(size(T,2),1); % sum over mode 2
    T = reshape(T, size_output);
elseif nargin == N+1
    U = cellfun(@(u,s) u(s,:), U, varargin, 'UniformOutput', false);
    T = cpdgen(U);
else
    error('cpdgen:index', ...
          'Either linear or subscripts indices should be provided');    
end
