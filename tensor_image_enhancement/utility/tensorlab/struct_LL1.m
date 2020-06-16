function [x,state] = struct_LL1(z,task,L)
%STRUCT_LL1 Structure of third factor matrix in a LL1 decomposition.
%   [x,state] = struct_LL1(z,[],L) generates x as
%
%      [z(:,ones(1,L(1))) ... z(:,length(L)*ones(1,L(end)))]
%
%   which is the structure of the third factor matrix of a rank-(Lr,Lr,1)
%   block term decomposition when written in CPD form. The structure state
%   stores information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_LL1(z,task,size_mat,pre,post) computes the right or left
%   Jacobian-vector product of this transformation, depending on the
%   structure task. Use the structure state and add the field 'r' of the
%   same shape as z or the field 'l' of the same shape as x to obtain the
%   structure task for computing the right and left Jacobian-vector
%   products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
if nargin < 3 || isempty(L), L = ones(1,size(z,2)); end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    R = sum(L);
    state.idx = sum(bsxfun(@gt,1:R,cumsum(L)'),1)+1;
    state.sum = sparse(1:R,state.idx,1);
    x = z(:,state.idx);
elseif ~isempty(task.r)
    x = task.r(:,task.idx);
    state = [];
elseif ~isempty(task.l)
    x = task.l*task.sum;
    state = [];
end

end
