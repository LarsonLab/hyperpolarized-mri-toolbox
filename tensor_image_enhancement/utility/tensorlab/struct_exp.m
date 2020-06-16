function [x,state] = struct_exp(z,task,t)
%STRUCT_EXP Matrix with columns as exponentials.
%   [x,state] = struct_exp(z,[],t) computes a matrix x in which the
%   jth column is equal to the exponential z(j)^t evaluated at the given
%   points t. z should be a row vector.
%
%   struct_exp(z,task,t) computes the right or left Jacobian-vector
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

%   Author:  Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2014/04/01   OD      Initial version

state = [];

if nargin < 2, task = []; end
if nargin < 3
    error('struct_exp:t','Please supply evaluation points.');
end
if numel(z) ~= length(z) || size(z,2) ~= length(z)
    error('struct_exp:z','z should be a row vector.');
end

if isempty(task)
    x = function_exp(z,t);
    state.derpoles = bsxfun(@times,function_exp(z,t-1),t(:));
elseif ~isempty(task.r)
    x = bsxfun(@times,task.derpoles,task.r);
elseif ~isempty(task.l)
    x = sum(conj(task.derpoles).*task.l,1);
end

end

function f = function_exp(poles,t)
    f = bsxfun(@power,poles,t(:));
end