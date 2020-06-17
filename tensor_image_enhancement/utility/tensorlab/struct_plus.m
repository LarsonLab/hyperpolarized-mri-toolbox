function [x,state] = struct_plus(z,task,cnst)
%STRUCT_PLUS Plus.
%   [x,state] = struct_plus(z,[],cnst) computes x as bsxfun(@plus,z,cnst). If z
%   and cnst are of the same size, this is equivalent to x = z+cnst. If z is a
%   cell, each of the matrices in the cell should have the same size, and x is
%   computed as bsxfun(@plus, sum(cat(3, z{:}), 3), cnst), i.e. the matrices in
%   z are summed and added to the constant cnst. The structure state stores
%   information which is reused in computing the right and left Jacobian-vector
%   products.
%
%   [x,state] = struct_plus(z, []) is equivalent to struct_plus(z, [], 0).
%
%   struct_plus(z,task,cnst) computes the right or left Jacobian-vector
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
%   See also struct_times.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
if nargin < 3, cnst = 0; end
state = [];

if iscell(z)
    N = ndims(z{1});
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        x = bsxfun(@plus,sum(cat(N+1,z{:}),N+1),cnst);
    elseif ~isempty(task.r)
        x = sum(cat(N+1, task.r{:}),N+1);
    elseif ~isempty(task.l)
        x = repmat({task.l}, 1, length(z));
    end
else % 
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        x = bsxfun(@plus,z,cnst);
    elseif ~isempty(task.r)
        x = task.r;
    elseif ~isempty(task.l)
        x = task.l;
    end
end
end
