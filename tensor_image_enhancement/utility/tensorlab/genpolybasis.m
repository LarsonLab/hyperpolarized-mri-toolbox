function y = genpolybasis(t,degree,type,barycentric)
%GENPOLYBASIS Polynomial basis.
%   Y = GENPOLYBASIS(t,degree,type) computes a basis of polynomials of
%   given type up to the given degree in the given function evaluation
%   points t. The following bases are supported:
%       - monomial          Monomial basis
%       - chebyshev         Chebyshev basis of the first kind
%       - chebyshev2        Chebyshev basis of the second kind
%       - legendre          Legendre basis
%   The columns are ranked with decreasing degree, analogous to the matlab
%   polyval script.
%
%   Author:  Otto Debals         (otto.debals@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version history:
% - 2014/03/31    OD        Initial version
% - 2016/02/08    OD        Renewed version with Chebyshev2 and barycentric
% - 2016/02/08    OD        Renewed version with matrix t as points

if nargin<3, type = 'Monomial'; end
if nargin<4, barycentric = false; end
if size(t,1)==1, t = t.'; end

d = degree;
type = lower(type);

if size(t,2)==1
    if ~barycentric
        switch type
            case 'monomial'
                y = ones(numel(t),d+1);
                for i = d:-1:1, y(:,i) = y(:,i+1).*t; end
            case 'legendre'
                y = [zeros(numel(t),d) ones(numel(t),1)];
                if d>0, y(:,end-1) = t; end
                for i = d-1:-1:1
                    n = d-i;
                    y(:,i) = ((2*n+1)*t.*y(:,i+1)-n*y(:,i+2))/(n+1);
                end
            case 'chebyshev'
                y = [zeros(numel(t),d) ones(numel(t),1)];
                if d>0, y(:,end-1) = t; end
                for i = d-1:-1:1
                    y(:,i) = 2*t.*y(:,i+1)-y(:,i+2);
                end
            case 'chebyshev2'
                y = [zeros(numel(t),d) ones(numel(t),1)];
                if d>0, y(:,end-1) = 2*t; end
                for i = d-1:-1:1
                    y(:,i) = 2*t.*y(:,i+1)-y(:,i+2);
                end
            otherwise
                error('genpolybasis:type','Unsupported type of polynomial basis');
        end
    else
        switch type
            case 'monomial'
                ta = linspace(-1,1,d+1);
                w = zeros(1,d+1);
                for i = 0:d, w(i+1) = (-1)^i*nchoosek(d,i); end
            case 'legendre'
                error('genpolybasis:legendrebary','Barycentric option for Legendre basis not yet supported!');
            case 'chebyshev'
                ta = cos((2*(d:-1:0)+1)*(pi/(2*d+2)));
                w = (-1).^(d:-1:0).*sin((2*(d:-1:0)+1)*(pi./(2*d+2)));
            case 'chebyshev2'
                ta = cos(pi*(d:-1:0)/d);
                w = (-1).^(0:d);
                w(1) = 1/2; w(end) = w(end)/2;
            otherwise
                error('genpolybasis:type','Unsupported type of polynomial basis');
        end
        P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
        y = bsxfun(@rdivide,P,sum(P,2));
        y(isnan(y)) = 1;
    end
else
    if ~barycentric
        switch type
            case 'monomial'
                y = ones([size(t) d+1]);
                for i = d:-1:1, y(:,:,i) = y(:,:,i+1).*t; end
            case 'legendre'
                y = ones([size(t) d+1]);
                if d>0, y(:,:,end-1) = t; end
                for i = d-1:-1:1
                    n = d-i;
                    y(:,:,i) = ((2*n+1)*t.*y(:,:,i+1)-n*y(:,:,i+2))/(n+1);
                end
            case 'chebyshev'
                y = ones([size(t) d+1]);
                if d>0, y(:,:,end-1) = t; end
                for i = d-1:-1:1
                    y(:,:,i) = 2*t.*y(:,:,i+1)-y(:,:,i+2);
                end
            case 'chebyshev2'
                y = ones([size(t) d+1]);
                if d>0, y(:,:,end-1) = 2*t; end
                for i = d-1:-1:1
                    y(:,:,i) = 2*t.*y(:,:,i+1)-y(:,:,i+2);
                end
            otherwise
                error('genpolybasis:type','Unsupported type of polynomial basis');
        end
    else
        switch type
            case 'monomial'
                ta = linspace(-1,1,d+1);
                w = zeros(1,d+1);
                for i = 0:d, w(i+1) = (-1)^i*nchoosek(d,i); end
            case 'legendre'
                error('genpolybasis:legendrebary','Barycentric option for Legendre basis not yet supported!');
            case 'chebyshev'
                ta = cos((2*(d:-1:0)+1)*(pi/(2*d+2)));
                w = (-1).^(d:-1:0).*sin((2*(d:-1:0)+1)*(pi./(2*d+2)));
            case 'chebyshev2'
                ta = cos(pi*(d:-1:0)/d);
                w = (-1).^(0:d);
                w(1) = 1/2; w(end) = w(end)/2;
            otherwise
                error('genpolybasis:type','Unsupported type of polynomial basis');
        end
        w = reshape(w,[1 1 numel(w)]);
        ta = reshape(ta,[1 1 numel(ta)]);
        P = bsxfun(@rdivide,w,bsxfun(@minus,t,ta));
        y = bsxfun(@rdivide,P,sum(P,3));
        y(isnan(y)) = 1;
    end
end

end