function T = ttgen(U,varargin)
%TTGEN Generates full tensor from TT format
%   T = ttgen(U) computes the tensor T as the product of two matrices and
%   length(U)-2 third order tensors (the tensor train format). 
%
%   T = ttgen(U,ind) computes only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = ttgen(U,i,j,k,...) computes only the elements T(i,j,k,...) from the
%   tensor. The number of indices should match the order of the tensor given by
%   length(U). The colon operator can be used as a string, e.g., for a third
%   order tensor, ttgen(U,i,':',k) can be used.
%
%   See also cpdgen, lmlragen, ll1gen, btdgen.
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/13   NV      Added incomplete ttgen 
% - 2015/12/17   NV      Initial version
    
    ttr = [1 cellfun('size', U(2:end), 1) 1];
    size_tens = cellfun('size', U, 2);
    size_tens(1) = size(U{1},1);
    N = length(size_tens);
    
    if nargin == 1
        T = U{1};
        for n = 2:length(U)
            T = T*reshape(U{n},ttr(n),size_tens(n)*ttr(n+1));
            T = reshape(T, prod(size_tens(1:n)), ttr(n+1));
        end
        T = reshape(T, size_tens);
    elseif nargin == 2
        if ischar(varargin{1}) && strcmp(varargin{1},':')
            T = reshape(ttgen(U), [], 1);
            return;
        end
        size_output = size(varargin{1});
        if nargin == 2 % linear indexing
            sub = cell(1, N);
            [sub{:}] = ind2sub(size_tens, varargin{1}(:));
        end
        U(2:end) = cellfun(@(u) permute(u, [2 1 3]), U(2:end), 'UniformOutput', false);
        T = U{1}(sub{1},:);
        for n = 2:N
            T = bsxfun(@times, T, U{n}(sub{n},:,:));
            T = squeeze(sum(T, 2));
        end
        T = reshape(T, size_output);
    elseif nargin == N+1
        U{1} = U{1}(varargin{1},:);
        U(2:end) = cellfun(@(u,s) u(:,s,:), U(2:end), varargin(2:end), ...
                           'UniformOutput', false);
        T = ttgen(U);
    else 
        error('ttgen:index', ...
              'Either linear or subscripts indices should be provided');    
    end
    
end
