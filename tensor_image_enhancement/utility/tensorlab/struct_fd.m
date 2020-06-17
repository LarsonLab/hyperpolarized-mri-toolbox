function [x,state] = struct_fd(z,task,order,t)
%STRUCT_FD Finite differences 
%   [x,state] = STRUCT_FD(z,[]) computes the first-order finite differeces for
%   each column of z using a step length of one. The structure state stores
%   information which is reused in computing the right and left Jacobian-vector
%   products.
%   
%   [x,state] = struct_fd(z,[],order) computes the finite differences of a
%   specified order. Currently, first-order (order=1) and second-order
%   (order=2) are supported.
%
%   [x,state] = struct_fd(z,[],order,t) computes the finites differences as
%   if z is evaluated in points t. The points t can be a column vector if the
%   points are equal for all columns of z, or a matrix of size(z) if each
%   column has its own evaulation points. If t is a scalar, the step length
%   is assumed to be constant and equal to 1/t. 
%
%   struct_fd(z,task,order,t) computes the right or left Jacobian-vector product
%   of this transformation, depending on the structure task. Use the structure
%   state and add the field 'r' of the same shape as z or the field 'l' of the
%   same shape as x to obtain the structure task for computing the right and
%   left Jacobian-vector products
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
%
% Version History:
% - 2015/01/23   NV      Initial version
    
    if nargin < 3, order = 1; end
    if nargin < 4, t = 1; end
    if numel(t) > 1
        if size(t, 1) == 1, t = t.'; end
        switch order
          case 1
            h = 1./diff(t, 1, 1);
          case 2
            h = 1./(diff(t(1:end-1,:),1,1).*diff(t(2:end,:),1,1));
        end 
    else 
        h = t.^(-order);
    end
    state = [];
    
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        x = bsxfun(@times, h, diff(z, order, 1));
    elseif ~isempty(task.r)
        x = bsxfun(@times, h, diff(task.r, order, 1));
    elseif ~isempty(task.l)
        R = size(task.l, 2);
        task.l = bsxfun(@times, h, task.l);
        switch order
          case 1
            x = [-task.l; zeros(1, R)] + [zeros(1, R); task.l];
          case 2
            x = [task.l; zeros(2, R)] + [zeros(1,R); -2*task.l; zeros(1,R)] + ...
                [zeros(2, R); task.l];
          otherwise
            error('struct_fd:notImplementend',...
                  'Finite differences of order > 2 are not implemented yet.');
        end 
    end
    
end
