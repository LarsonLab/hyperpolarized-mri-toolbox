function [x,state] = struct_times(z,task,cnst)
%STRUCT_TIMES Times.
%   [x,state] = struct_times(z,[],cnst) computes x as
%   bsxfun(@times,z,cnst). If z and cnst are of the same size, this is
%   equivalent to x = z.*cnst. The structure state stores information which
%   is reused in computing the right and left Jacobian-vector products.
%
%   struct_times(z,task,cnst) computes the right or left Jacobian-vector
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
%   
%   See also struct_plus.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2016/03/16   NV      Extended to cell of variables (matrices/tensors) +
%                        constant
% - 2014/02/01   LS      Initial version for matrix and constant
    

if nargin < 2, task = []; end
if nargin < 3, cnst = 1; end
state = [];

right = isstruct(task) && isfield(task,'r') && ~isempty(task.r);
left  = isstruct(task) && isfield(task,'l') && ~isempty(task.l);

if iscell(z)
    N = ndims(z{1});
    if ~left && ~right
        tmp = cat(N+1, z{:});
        x = bsxfun(@times,prod(tmp,N+1),cnst);
    elseif right
        tmp = cat(N+1, z{:});
        tmp = reshape(tmp, [], length(z));
        x = zeros(size(tmp,1),1);
        for n = 1:length(z)
            x = x + prod(tmp(:,[1:n-1 n+1:end]),2).*task.r{n}(:);
        end
        x = reshape(x, size(z{1}));
        x = bsxfun(@times, x, cnst);
    elseif left
        tmp = cat(N+1, z{:});
        tmp = conj(reshape(tmp, [], length(z)));
        x = cell(1, length(z));
        for n = 1:length(z)
            x{n} = prod(tmp(:,[1:n-1 n+1:end]),2).*task.l(:);
            x{n} = reshape(x{n}, size(z{1}));
            x{n} = bsxfun(@times, x{n}, conj(cnst));
        end
    end
else 
    if ~left && ~right
        x = bsxfun(@times,z,cnst);
    elseif right
        x = bsxfun(@times,task.r,cnst);
    elseif left
        x = bsxfun(@times,task.l,conj(cnst));
    end
end
