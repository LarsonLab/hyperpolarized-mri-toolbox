function [x,state] = struct_cauchy(z,task,u,v)
%STRUCT_CAUCHY Cauchy matrix.
%   [x,state] = STRUCT_CAUCHY(z) uses the generator vectors z{1} and z{2}
%   to compute a Cauchy matrix x with
%
%      x(i,j) = 1/(z{1}(i) - z{2}(j))
%
%   The structure state stores information which is reused in computing the
%   right and left Jacobian-vector products.
%
%   [x,state] = STRUCT_CAUCHY(z,[],u) uses the generator vector z and fixed
%   u to compute a Cauchy matrix x with
%
%      x(i,j) = 1/(u(i) - z(j))
%
%   [x,state] = STRUCT_CAUCHY(z,[],[],v) uses the generator vector z and
%   fixed v to compute a Cauchy matrix x with
%
%      x(i,j) = 1/(z(i) - v(j))
%
%   STRUCT_CAUCHY(z,task,u,v) computes the right or left Jacobian-vector
%   product of this transformation, depending on the structure task. Use
%   the structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and (dF(:)/dz(:).')'*task.l(:) +
%      conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.
%   
%   See also struct_hankel, struct_toeplitz.

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

if nargin < 2, task = []; end
if ~iscell(z), z = {z}; wascell = false; else wascell = true; end
if nargin<3, u = []; end
if nargin<4, v = []; end

if numel(z)>2, error('struct_cauchy:z','z should contain at most two variables');
elseif numel(z)==2 && nargin>2
    error('struct_cauchy:uv','If z contains two variables, u and v should not be passed.');
elseif numel(z)==1
    if nargin>3 && ~isempty(u)
        error('struct_cauchy:zuv','If z contains one variable, only one of u and v should be passed.');
    end
end

if ~isempty(u) && isrow(u), u = u.'; end
if ~isempty(v) && isrow(v), v = v.'; end
if ~isempty(u) && ~iscolumn(u), error('struct_cauchy:usize','u should contain a vector'); end
if ~isempty(v) && ~iscolumn(v), error('struct_cauchy:vsize','v should contain a vector'); end
if ~all(cellfun(@isvector,z)), error('struct_cauchy:zmat','z should contain vectors'); end

if numel(z)==1
    if size(z{1},2)>1, z{1} = z{1}.'; end
    % error('struct_cauchy:zsize','z should contain column vectors');
    if isempty(v) && ~isempty(u)
        z{2} = z{1};
        z{1} = u;
    elseif ~isempty(v) && isempty(u)
        z{2} = v;
    else
        error('struct_cauchy:z','z should contain two variables, or u/v should be given (not both)');
    end
else
    if isrow(z{1}), z{1} = z{1}.'; end
    % error('struct_cauchy:zsize','z should contain column vectors'); end
    if isrow(z{2}), z{2} = z{2}.'; end
end

state = [];
if isempty(task) || (isempty(task.l) && isempty(task.r))
    x = 1./bsxfun(@minus,z{1},z{2}.');
    state.x2 = x.^2;
elseif ~isempty(task.r)
    if iscell(task.r), taskr1 = task.r{1}; taskr2 = task.r{2};
    else taskr1 = task.r; taskr2 = task.r;
    end
    if ~isempty(u)
        x = bsxfun(@times,(task.x2),taskr2.');
    elseif ~isempty(v)
        x = bsxfun(@times,-(task.x2),taskr1);
    else
        x = bsxfun(@times,-(task.x2),taskr1)+bsxfun(@times,(task.x2),taskr2.');
    end
elseif ~isempty(task.l)
    tmp = conj(task.x2).*task.l;
    if ~isempty(u)
        x = sum(tmp,1).';
        if wascell, x = {x}; end
    elseif ~isempty(v)
        x = -sum(tmp,2);
        if wascell, x = {x}; end
    else 
        x{1} = -sum(tmp,2);
        x{2} = sum(tmp,1).';
    end
end

end