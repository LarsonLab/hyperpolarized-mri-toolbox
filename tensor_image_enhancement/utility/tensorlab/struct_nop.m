function [x,state] = struct_nop(z,task,varargin)
%STRUCT_NOP No operation. 
%   [x,state] = struct_nop(z,[],varargin) does not transform the variables. It
%   accepts a variable number of extra arguments, which are all passed on to
%   state.varargs which is a cell with one entry for each extra input argument.
%   This function is mainly useful for testing purposes.
%
%   struct_nop(z,task) computes the right or left Jacobian-vector
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

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

    state = struct;
    state.varargs = varargin;
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        x = z;
    elseif ~isempty(task.r)
        x = task.r;
    elseif ~isempty(task.l)
        x = task.l;
    end
end
