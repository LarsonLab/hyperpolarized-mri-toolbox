function [x,state] = struct_prod(z,task,dim)
%STRUCT_PROD Product of elements.
%   [x,state] = struct_prod(z,[],dim) computes x as squeeze(prod(z,dim)), or as
%   prod(z(:)) if dim is the empty matrix [] or not supplied.
%
%   struct_prod(z,task,dim) computes the right or left Jacobian-vector
%   product of this transformation, depending on the structure task. Use
%   the structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2016/03/16   NV      Cleaned code, improved performance and made it work
%                        for arbitrary order variables and dims. Removed
%                        constant factor A (use struct_times for that)                            
% - 2014/03/05   OD      Initial version

if nargin < 2, task = []; end
if nargin < 3, dim = []; end
state = [];

if ~isempty(dim) 
    if ~isscalar(dim)
        error('struct_prod:dim', 'dim should be a scalar');
    end
    if dim < 1 || dim > ndims(z)
        error('struct_prod:dim', 'dim should be >= 1 and <= ndims(z)');
    end
end

right = isstruct(task) && isfield(task,'r') && ~isempty(task.r);
left  = isstruct(task) && isfield(task,'l') && ~isempty(task.l);

% if dim is empty, treat as row vector
nodim = isempty(dim);
if nodim
    z = z(:).';
    dim = 2;
end

% compute sizes
sz = size(z);
odim = 1:ndims(z); % other dim
odim(dim) = [];
m = prod(sz(odim));
n = sz(dim);

if ~left && ~right
    if isempty(dim), x = prod(z(:)); 
    else x = squeeze(prod(z, dim)); end 
    % Precompute values for task.r and task.l
    state = struct;
    % permute and reshape to get a matrix with dim as column mode
    z = permute(z, [odim, dim]);
    z = reshape(z, m, n);
    % expand z to get product tensor
    state.reshapedz = z;
    state.expandedz = z(:,:,ones(1,n));
    % construct indicies of diagonal mode-1 vectors
    idx = linspace(1,n^2,n)*m - m;
    state.diagidx = bsxfun(@plus, (1:m).', idx);
elseif right
    % permute and reshape to get a matrix with dim as column mode
    if nodim
        task.r = task.r(:).';
    else 
        task.r = permute(task.r, [odim, dim]);
        task.r = reshape(task.r, m, n);
    end
    % expand z to get product tensor
    tmp = task.expandedz;
    idx = task.diagidx;
    % replace z by task.r in the right places
    tmp(idx) = tmp(idx) - (task.reshapedz - task.r);
    % compute derivative
    x = sum(prod(tmp, 2), 3);
    % fix dimensions
    if length(odim) > 1, x = reshape(x, sz(odim));
    elseif dim == 1, x = x.'; end
elseif left
    % extract precomputed values
    tmp = task.expandedz; 
    idx = task.diagidx;
    % replace z by task.l in the right places
    tmp(idx) = tmp(idx) - bsxfun(@minus, task.reshapedz, conj(task.l(:)));
    % compute derivative
    x = conj(prod(tmp, 2));
    % fix dimensions
    x = reshape(x, [sz(odim), sz(dim)]);
    x = permute(x, [1:dim-1 length(sz) dim:length(sz)-1]);
end

end
