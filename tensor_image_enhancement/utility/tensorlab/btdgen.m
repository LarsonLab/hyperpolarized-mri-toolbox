function T = btdgen(U,varargin)
%BTDGEN Generate full tensor given a BTD.
%   T = btdgen(U) computes the tensor T as the sum of the R block terms
%   U{r} which are computed as tmprod(U{r}{N+1},U{r}(1:N),1:N), where N is
%   the order of the tensor T.
%
%   T = btdgen(U,ind) computes only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = btdgen(U,i,j,k,...) computes only the elements T(i,j,k,...) from the
%   tensor. The number of indices should match the order of the tensor given
%   by length(U{1}). The colon operator can be used as a string, e.g., for a
%   third order tensor, btdgen(U,i,':',k) can be used. 
%
%   See also cpdgen, lmlragen, ll1gen, ttgen, btdres.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/01/13  NV   Added incomplete btdgen

N = length(U{1})-1;
size_tens = cellfun('size',U{1}(1:N),1);
if nargin == 1
    T = reshape(U{1}{1}*mtkronprod(U{1}{end},U{1}(1:N),1,'H'),size_tens);
    for r = 2:length(U)
        T = T+reshape(U{r}{1}*mtkronprod(U{r}{end},U{r}(1:N),1,'H'),size_tens);
    end
elseif nargin == 2
    if ischar(varargin{1}) && strcmp(varargin{1},':')
        T = reshape(btdgen(U), [], 1);
        return;
    end
    size_output = size(varargin{1});
    if nargin == 2 % linear indexing
        sub = cell(1, N);
        [sub{:}] = ind2sub(size_tens, varargin{1}(:));
    end
    T = 0;
    for r = 1:length(U)
        S = U{r}{N+1};
        size_core = cellfun('size',U{r}(1:end-1),2);
        idx = cell(1,length(U{r})-1);
        for i = 1:numel(S)
            [idx{:}] = ind2sub(size_core,i);
            tmp = S(idx{:})*U{r}{1}(sub{1},idx{1});
            for n = 2:length(size_core)
                tmp = tmp .*U {r}{n}(sub{n},idx{n});
            end
            T = T + tmp;
        end
    end
    T = reshape(T, size_output);
elseif nargin == N + 1
    for r = 1:length(U)
        U{r}(1:N) = cellfun(@(u,s) u(s,:), U{r}(1:N), varargin, 'UniformOutput', ...
                            false);
    end
    T = btdgen(U);
else 
    error('btdgen:index', ...
          'Either linear or subscripts indices should be provided');    
end
