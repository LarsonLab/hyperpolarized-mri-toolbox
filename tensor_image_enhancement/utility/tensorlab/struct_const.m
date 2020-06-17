function [x,state] = struct_const(z,task,mask)
%STRUCT_CONST Keep parts of z constant.
%   [x,state] = struct_const(z,[]) keeps the factor z constant by setting
%   derivatives to zero. Note that the underlying variables can still be
%   updated by factors. 
%
%   [x,state] = struct_const(z,[],mask) keeps the entries for which the binary
%   mask == 0 constant. mask is a binary matrix of size size(z). 
%
%   struct_const(z,task,mask) computes the right or left Jacobian-vector
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
    
%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
    
    if nargin < 2, task = []; end
    if nargin < 3, 
        mask = zeros(size(z)); 
    elseif any(size(z) ~= size(mask))
        error('struct_const:dimensionMismatch', ...
              'size(mask) should match size(z)');
    end
    
    state = [];
    
    if isempty(task) || (isempty(task.r) && isempty(task.l))
        x = z;
    elseif ~isempty(task.r)
        x = mask.*task.r;
    elseif ~isempty(task.l)
        x = mask.*task.l;
    end
end
