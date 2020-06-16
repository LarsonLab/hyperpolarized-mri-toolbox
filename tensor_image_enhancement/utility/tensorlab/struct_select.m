function [x,state] = struct_select(z,task, n)
%STRUCT_SELECT Select entry from cell variable z.
%   [x,state] = struct_select(z, [], n) selects the nth matrix from a cell of
%   variables z. The structure state stores information which is reused in
%   computing the right and left Jacobian-vector products.
%
%   struct_select(z,task,n) computes the right or left Jacobian-vector product
%   of this transformation, depending on the structure task. Use the structure
%   state and add the field 'r' of the same shape as z or the field 'l' of the
%   same shape as x to obtain the structure task for computing the right and
%   left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial derivative
%   which treats conj(z) (z) as constant. The output has the same shape as x or
%   z for the right and left Jacobian-vector products, respectively.
%   
%   See also struct_times.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
    
if nargin < 2, task = []; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    x = z{n};
elseif ~isempty(task.r)
    x = task.r{n};
elseif ~isempty(task.l)
    x = cellfun(@(u) zeros(size(u)), z, 'UniformOutput', false);
    x{n} = task.l;
end

end