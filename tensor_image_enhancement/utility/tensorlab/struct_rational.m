function [x,state] = struct_rational(z,task,t,basis,barycentric)
%STRUCT_RATIONAL Matrix with columns as rational functions.
%   [x,state] = STRUCT_RATIONAL(z,[],t) computes a matrix x in which the
%   jth column is equal to the rational function
%
%      polyval(z{1}(j,:),s)./polyval([z{2}(j,:) 1],s)
%
%   evaluated at the points s, defined as
%
%      (t-0.5*(min(t)+max(t)))/(0.5*(max(t)-min(t))).
%
%   The degree of the numerator and denomator are equal to size(z{1},2)-1
%   and size(z{2},2), respectively. The structure state stores information
%   which is reused in computing the right and left Jacobian-vector
%   products.
%
%   [x,state] = STRUCT_RATIONAL(z,[],t,basis,barycentric) uses the given
%   basis, which can be one of the following:
%       - monomial          Monomial basis (default)
%       - chebyshev         Chebyshev basis of the first kind
%       - chebyshev2        Chebyshev basis of the second kind
%       - legendre          Legendre basis
%   If barycentric is true, a barycentric interpolation is used (default:
%   false). Instead of using the coefficients of the corresponding basis,
%   function values on a specific grid are used, depending of the chosen
%   basis. This should improve numerical properties.
%
%   [x,state] = STRUCT_RATIONAL(z,[]) computes a matrix x in which the jth
%   column is equal to the rational function
%
%       polyval(z{1}(j,:),z{3})./polyval([z{2}(j,:) 1],z{3})
%
%   Hence, instead of using fixed points t, the points are optimized as
%   well.
%
%   STRUCT_RATIONAL(z,task,t) computes the right or left Jacobian-vector
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
%   See also struct_poly, struct_rbf.

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
if ~iscell(z),
    error('struct_rational:z','z should be a cell with two or three elements.');
end
if nargin < 3 && numel(z) < 3
    error('struct_rational:t','Please supply evaluation points.');
end
if numel(z) > 2 && nargin > 2 && ~isempty(t)
    error('struct_rational:zt',...
        'The points t should not be given twice, both explicitly and in z');
end
if nargin < 3 || (nargin > 2 && isempty(t))
    optimizepoints = true;  % Optimize the point set t
    t = z{3};
else optimizepoints = false;
end
if nargin < 4, basis = 'monomial'; end
if nargin < 5, barycentric = false; end
if optimizepoints && barycentric,
    error('struct_rational:optimizeandbary','Not yet supported!');
end
if size(t,1)==1, t = t.'; end

a = z{1};
b = z{2};

state = [];

if ~isstruct(task) || ~isfield(task,'persistent')
    switch basis
        case 'monomial'
            if barycentric, basisf = @monomial_bary_basis;
            else basisf = @monomial_basis;
            end
        case 'chebyshev'
            if barycentric, basisf = @chebyshev_bary_basis;
            else basisf = @chebyshev_basis;
            end
        case 'chebyshev2'
            if barycentric, basisf = @chebyshev2_bary_basis;
            else basisf = @chebyshev2_basis;
            end
        case 'legendre'
            if barycentric, error('struct_rational:legendrebary',...
                    'Barycentric interpolation for Legendre basis not yet supported!'); end
            basisf = @legendre_basis;
        otherwise
            error('struct_rational:basis','Not yet supported!')
    end
end

if ~optimizepoints
    % Not optimizing the points t
    
    if ~isstruct(task) || ~isfield(task,'persistent')
        % The basis is created only in the first call
        mx = max(t,[],1);
        mn = min(t,[],1);
        t = bsxfun(@rdivide,bsxfun(@minus,t,.5*(mn+mx)),.5*(mx-mn));
        t(:,mn==mx) = ones(size(t,1),sum(mn==mx));
        
        if ~barycentric
            if size(a,2)>=size(b,2)+1
                Ba = basisf(t,size(a,2)-1);
                if size(t,2)==1, Bb = Ba(:,end-size(b,2):end);
                else Bb = Ba(:,:,end-size(b,2):end); end
            else
                Bb = basisf(t,size(b,2));
                if size(t,2)==1, Ba = Bb(:,end-size(a,2)+1:end);
                else Ba = Bb(:,:,end-size(a,2)+1:end); end
            end
        else
            Ba = basisf(t,size(a,2)-1);
            Bb = basisf(t,size(b,2));
        end
        state.persistent = {Ba,Bb};
    else
        Ba = task.persistent{1};
        Bb = task.persistent{2};
    end
    
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        x = eval(Ba,a);
        if ~isempty(b)
            state.denom = eval(Bb,[b ones(size(b,1),1)]);
            x = x./state.denom;
            state.deriv = -x./state.denom;
        end
    elseif ~isempty(task.r)
        x = eval(Ba,task.r{1});
        if ~isempty(b)
            c = [task.r{2} zeros(size(b,1),1)];
            x = x./task.denom + ...
                task.deriv.*eval(Bb,c);
        end
    elseif ~isempty(task.l)
        x = {zeros(size(a)),zeros(size(b))};
        if size(t,2)==1, Baext = reshape(Ba,[size(Ba,1) 1 size(Ba,2)]);
        else Baext = Ba;
        end
        
        if isempty(b)
            x{1} = squeeze(sum(bsxfun(@times,conj(Baext),task.l),1));
        else
            if size(t,2)==1, Bbext = reshape(Bb,[size(Bb,1) 1 size(Bb,2)]);
            else Bbext = Bb;
            end
            x{1} = squeeze(sum(bsxfun(@rdivide,conj(Baext),conj(task.denom)./task.l),1));
            Bcext = Bbext(:,:,1:end-1);
            x{2} = squeeze(sum(bsxfun(@times,conj(Bcext),conj(task.deriv).*task.l),1));
        end
    end
else
    % Optimizing the points t
    
    if isempty(task) || (isempty(task.l) && isempty(task.r))
        if size(a,2)>=size(b,2)+1
            [Ba,Bad] = basisf(t,size(a,2)-1);
            if size(t,2)==1, Bb = Ba(:,end-size(b,2):end);
                Bbd = Bad(:,end-size(b,2):end);
            else Bb = Ba(:,:,end-size(b,2):end);
                Bbd = Bad(:,:,end-size(b,2):end);
            end
        else
            [Bb,Bbd] = basisf(t,size(b,2));
            if size(t,2)==1, Ba = Bb(:,end-size(a,2)+1:end);
                Bad = Bbd(:,end-size(a,2)+1:end);
            else Ba = Bb(:,:,end-size(a,2)+1:end);
                Bad = Bbd(:,:,end-size(a,2)+1:end);
            end
        end
        state.Ba = Ba;
        state.Bb = Bb;
        
        x = eval(Ba,a);
        state.dera = eval(Bad,a);
        if ~isempty(b)
            state.denom = eval(Bb,[b ones(size(b,1),1)]);
            x = x./state.denom;
            state.deriv = -x./state.denom;
            state.derb = eval(Bbd,[b ones(size(b,1),1)]);
        end
        
    elseif ~isempty(task.r)
        Ba = task.Ba;
        x = eval(Ba,task.r{1});
        if ~isempty(b)
            c = [task.r{2} zeros(size(b,1),1)];
            x = x./task.denom + ...
                task.deriv.*eval(task.Bb,c);
            x = x + bsxfun(@times,task.r{3},task.dera./task.denom + task.deriv.*task.derb);
        else
            x = x + bsxfun(@times,task.dera,task.r{3});
        end
        
    elseif ~isempty(task.l)
        Ba = task.Ba;
        x = {zeros(size(z{1})),zeros(size(z{2})),zeros(size(t))};
        if size(t,2)==1, Baext = reshape(Ba,[size(Ba,1) 1 size(Ba,2)]);
        else Baext = Ba;
        end
        if isempty(b)
            x{1} = squeeze(sum(bsxfun(@times,conj(Baext),task.l),1));
            x{3} = conj(task.dera).*task.l;
        else
            Bb = task.Bb;
            if size(t,2)==1, Bbext = reshape(Bb,[size(Bb,1) 1 size(Bb,2)]);
            else Bbext = Bb;
            end
            x{1} = squeeze(sum(bsxfun(@rdivide,conj(Baext),conj(task.denom)./task.l),1));
            Bcext = Bbext(:,:,1:end-1);
            x{2} = squeeze(sum(bsxfun(@times,conj(Bcext),conj(task.deriv).*task.l),1));
            x{3} = conj(task.dera./task.denom+task.deriv.*task.derb).*task.l;
        end
        if size(t,2)==1, x{3} = sum(x{3},2); end
    end
end

end

function f = eval(B,c)
if ~ismatrix(B)
    f = sum(bsxfun(@times,B,reshape(c,[1 size(c)])),3);
else
    f = B*c.';
end
end

function [T,Td] = monomial_basis(t,d)
% Setup of monomial basis
if size(t,2)==1
    T = bsxfun(@power,t,d:-1:0);
    if nargout>1
        Td = zeros(size(T));
        Td(:,1:end-2) = bsxfun(@times,d:-1:2,T(:,2:end-1));
        Td(:,end-1) = ones(size(Td,1),1);
    end
else
    T = bsxfun(@power,t,reshape(d:-1:0,[1 1 d+1]));
    if nargout>1
        Td = zeros(size(T));
        Td(:,:,1:end-2) = bsxfun(@times,reshape(d:-1:2,[1 1 d-1]),T(:,:,2:end-1));
        Td(:,:,end-1) = ones(size(t));
    end
end
end

function T = monomial_bary_basis(t,d)
% Setup of monomial basis using monomial formula
if size(t,2)==1
    ta = linspace(-1,1,d+1);
    w = zeros(1,d+1);
    for i = 0:d, w(i+1) = (-1)^i*nchoosek(d,i); end
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,2));
    T(isnan(T)) = 1;
else
    ta = reshape(linspace(-1,1,d+1),[1 1 d+1]);
    w = zeros(1,1,d+1);
    for i = 0:d, w(i+1) = (-1)^i*nchoosek(d,i); end
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,3));
    T(isnan(T)) = 1;
end
end

function [T,Td] = chebyshev_basis(t,d)
% Setup of Chebyshev basis of the first kind

% Possible to use explicit function, but then problems for t in [-Inf,-1]
% and [1,Inf] because of 1/sqrt(1-t^2)?
% T = cos(acos(t(:))*(d:-1:0));
if size(t,2)==1
    T = [zeros(numel(t),d) ones(size(t,1),1)];
    Td = zeros(numel(t),d+1);
    if d>0, T(:,end-1) = t(:); Td(:,end-1) = ones(size(t(:))); end
    for i = d-1:-1:1
        T(:,i) = 2*t(:).*T(:,i+1)-T(:,i+2);
        if nargout>1, Td(:,i) = 2*T(:,i+1) + 2.*t(:).*Td(:,i+1)-Td(:,i+2); end
    end
else
    T = ones([size(t) d+1]);
    Td = zeros([size(t) d+1]);
    if d>0, T(:,:,end-1) = t; Td(:,:,end-1) = ones(size(t)); end
    for i = d-1:-1:1
        T(:,:,i) = 2*t.*T(:,:,i+1)-T(:,:,i+2);
        if nargout>1, Td(:,:,i) = 2*T(:,:,i+1) + 2.*t.*Td(:,:,i+1)-Td(:,:,i+2); end
    end
end
end

function T = chebyshev_bary_basis(t,d)
% Setup of Chebyshev basis of the first kind using barycentric formula
if size(t,2)==1
    ta = cos((2*(d:-1:0)+1)*(pi/(2*d+2)));
    w = (-1).^(d:-1:0).*sin((2*(d:-1:0)+1)*(pi./(2*d+2)));
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,2));
    T(isnan(T)) = 1;
else
    ta = reshape(cos((2*(d:-1:0)+1)*(pi/(2*d+2))),[1 1 d+1]);
    w = (-1).^(d:-1:0).*sin((2*(d:-1:0)+1)*(pi./(2*d+2)));
    w = reshape(w,[1 1 d+1]);
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,3));
    T(isnan(T)) = 1;
end
end

function [T,Td] = chebyshev2_basis(t,d)
% Setup of Chebyshev basis of the second kind

% Possible to use explicit function, but then problems for t in [-Inf,-1]
% and [1,Inf] because of 1/sqrt(1-t^2)?
% T = sin(acos(t(:))*(d+1:-1:1))./sin(acos(t(:)));
if size(t,2)==1
    T = [zeros(numel(t),d) ones(size(t,1),1)];
    Td = zeros(numel(t),d+1);
    if d>0, T(:,end-1) = 2*t(:); Td(:,end-1) = 2*ones(size(t(:))); end
    for i = d-1:-1:1
        T(:,i) = 2*t(:).*T(:,i+1)-T(:,i+2);
        if nargout>1, Td(:,i) = 2*T(:,i+1) + 2.*t(:).*Td(:,i+1)-Td(:,i+2); end
    end
else
    T = ones([size(t) d+1]);
    Td = zeros([size(t) d+1]);
    if d>0, T(:,:,end-1) = 2*t; Td(:,:,end-1) = 2*ones(size(t)); end
    for i = d-1:-1:1
        T(:,:,i) = 2*t.*T(:,:,i+1)-T(:,:,i+2);
        if nargout>1, Td(:,:,i) = 2*T(:,:,i+1) + 2.*t.*Td(:,:,i+1)-Td(:,:,i+2); end
    end
end
end

function T = chebyshev2_bary_basis(t,d)
% Setup of Chebyshev basis of the second kind using barycentric formula
if size(t,2)==1
    ta = cos(pi*(d:-1:0)/d);
    w = (-1).^(0:d);
    w(1) = 1/2; w(end) = w(end)/2;
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,2));
    T(isnan(T)) = 1;
else
    ta = cos(pi*(d:-1:0)/d);
    w = (-1).^(0:d);
    ta = reshape(ta,[1 1 d+1]); w = reshape(w,[1 1 d+1]);
    w(1) = 1/2; w(end) = w(end)/2;
    P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
    T = bsxfun(@rdivide,P,sum(P,3));
    T(isnan(T)) = 1;
end
end

function [T,Td] = legendre_basis(t,d)
% Setup of Legendre basis
if size(t,2)==1
    T = [zeros(numel(t),d) ones(size(t,1),1)];
    Td = zeros(numel(t),d+1);
    if d>0, T(:,end-1) = t(:); Td(:,end-1) = ones(size(t(:))); end
    for i = d-1:-1:1
        n = d-i;
        T(:,i) = ((2*n+1)*t(:).*T(:,i+1)-n*T(:,i+2))/(n+1);
        if nargout>1, Td(:,i) = ((2*n+1)*T(:,i+1) + (2*n+1).*t(:).*Td(:,i+1)-n*Td(:,i+2))/(n+1); end
    end
else
    T = ones([size(t) d+1]);
    Td = zeros([size(t) d+1]);
    if d>0, T(:,:,end-1) = t; Td(:,:,end-1) = ones(size(t)); end
    for i = d-1:-1:1
        n = d-i;
        T(:,:,i) = ((2*n+1)*t.*T(:,:,i+1)-n*T(:,:,i+2))/(n+1);
        if nargout>1, Td(:,:,i) = ((2*n+1)*T(:,:,i+1) + (2*n+1).*t.*Td(:,:,i+1)-n*Td(:,:,i+2))/(n+1); end
    end
end
end
