function [x,state] = struct_poly(z,task,t,basis,barycentric)
%STRUCT_POLY Matrix with columns as polynomials.
%   [x,state] = STRUCT_POLY(z,[],t) computes a matrix x in which the
%   jth column is equal to the polynomial
%
%      polyval(z(j,:),s)
%
%   evaluated at the points s, defined as
%
%      (t-0.5*(min(t)+max(t)))/(0.5*(max(t)-min(t))).
%
%   The degree of the polynomial is equal to size(z,2)-1. The structure
%   state stores information which is reused in computing the right and
%   left Jacobian-vector products.
%
%   [x,state] = STRUCT_POLY(z,[],t,basis,barycentric) uses the given basis,
%   which can be one of the following:
%       - monomial          Monomial basis (default)
%       - chebyshev         Chebyshev basis of the first kind
%       - chebyshev2        Chebyshev basis of the second kind
%       - legendre          Legendre basis
%   If barycentric is true, a barycentric interpolation is used (default:
%   false). Instead of using the coefficients of the corresponding basis,
%   function values on a specific grid are used, depending of the chosen
%   basis. This should improve numerical properties.
%
%   [x,state] = STRUCT_POLY(z,[]) computes a matrix x in which the jth
%   column is equal to the polynomial
%
%       polyval(z{1}(j,:),z{2}).
%
%   Hence, instead of using fixed points t, the points are optimized as
%   well.
%
%   STRUCT_POLY(z,task,t) computes the right or left Jacobian-vector
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
%   See also struct_rational, struct_rbf.

%   Authors: Otto Debals (Otto.Debals@esat,kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%   [2] J.P. Berrut, L.N. Trefethen, "Barycentric Lagrange Interpolation,"
%       SIAM Review, Vol. 46, No. 3, pp. 501-517, 2004.

if nargin < 2, task = []; end
if nargin < 3 && (~iscell(z) || numel(z)<2)
    error('struct_poly:t','Please supply evaluation points.');
end
if (iscell(z) && numel(z) > 1) && (nargin > 2 && ~isempty(t))
    error('struct_poly:zt','The points t should not be given twice, both explicitly and in z');
end
if nargin < 4, basis = 'monomial'; end
if nargin < 5, barycentric = false; end

if nargin < 3 || isempty(t) % Optimize points t
    if isempty(task) || (isempty(task.r) && isempty(task.l))
        [x,state] = struct_rational({z{1},[],z{2}},task,[],basis,barycentric);
    elseif ~isempty(task.r)
        task.r = {task.r{1},[],task.r{2}};
        [x,state] = struct_rational({z{1},[],z{2}},task,[],basis,barycentric);
    elseif ~isempty(task.l)
        [x,state] = struct_rational({z{1},[],z{2}},task,[],basis,barycentric); x = x([1 3]);
    end
else
    if isempty(task) || (isempty(task.r) && isempty(task.l))
        [x,state] = struct_rational({z,[]},task,t,basis,barycentric);
    elseif ~isempty(task.r)
        task.r = {task.r,[]};
        [x,state] = struct_rational({z,[]},task,t,basis,barycentric);
    elseif ~isempty(task.l)
        [x,state] = struct_rational({z,[]},task,t,basis,barycentric); x = x{1};
    end
end

end
