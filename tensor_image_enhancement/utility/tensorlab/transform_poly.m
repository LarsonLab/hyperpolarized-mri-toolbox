function b = transform_poly(a,polya,polyb)
%TRANSFORM_POLY Transform a polynomial
%   coeff = transform_poly(a,polya,polyb) transforms the coefficients a
%   (which can be a matrix consisting of multiple polynomials in each row)
%   from polya into coefficients b from polyb. Descending powers are used.
%   Each pola and polyb are structures with the following fields:
%       poly.type       - Type of polynomial, e.g. Monomial, Legendre,
%                         Chebyshev, Chebyshev2
%       poly.domain     - Domain of the polynomial coefficients, consisting
%                         of two elements, e.g. [0,1]

% Authors: Otto Debals         (Otto.Debals@esat.kuleuven.be)
%          Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/02/08    OD    Chebyshev basis of the second kind
% - 2014/04/01    OD    Multiple polynomials possible
% - 2014/03/31    OD    Initial version

if size(a,1) > 1
    a = mat2cell(a,ones(size(a,1),1),size(a,2));
    b = cell2mat(cellfun(@(x) transform_poly(x,polya,polyb),a,...
        'UniformOutput',false));
    return
end

if nargin < 3, polyb = struct('type','monomial','domain',[-1 1]); end

if any(polya.domain(:)~=polyb.domain(:))
    if ~strcmp(polya.type,polyb.type)
        if any(strcmp(polya.type,{'monomial','Monomial'}))
            %                 | Domain A | Domain B   |
            %  ---------------------------------------|
            %   Type A (mono) |    ------------>|     |
            %   -  -  -  -  - | -  -  -  | -  - |-  - |
            %   Type B        |          |      v     |
            %   -  -  -  -  - | -  -  -  | -  -  -  - |
            
            tmpstruct = struct('type',polya.type,'domain',polyb.domain);
            tmpcoeff = transform_poly(a,polya,tmpstruct);
            b = transform_poly(tmpcoeff,tmpstruct,polyb);
        elseif any(strcmp(polyb.type,{'monomial','Monomial'}))
            %                 | Domain A | Domain B   |
            %  ---------------------------------------|
            %   Type A        |    |     |            |
            %   -  -  -  -  - | -  |  -  | -  -  -  - |
            %   Type B (mono) |    v------------>     |
            %   -  -  -  -  - | -  -  -  | -  -  -  - |
            
            tmpstruct = struct('type',polyb.type,'domain',polya.domain);
            tmpcoeff = transform_poly(a,polya,tmpstruct);
            b = transform_poly(tmpcoeff,tmpstruct,polyb);
        else
            %                 | Domain A | Domain B   |
            %  ---------------------------------------|
            %   Type A        |    |     |            |
            %   -  -  -  -  - | -  |  -  | -  -  -  - |
            %   Type C (mono) |    v----------->|     |
            %   -  -  -  -  - | -  -  -  | -  - |-  - |
            %   Type B        |          |      v     |
            %   -  -  -  -  - | -  -  -  | -  -  -  - |
            
            tmpstruct1 = struct('type','monomial','domain',polya.domain);
            tmpcoeff1 = transform_poly(a,polya,tmpstruct1);
            tmpstruct2 = struct('type','monomial','domain',polyb.domain);
            tmpcoeff2 = transform_poly(tmpcoeff1,tmpstruct1,tmpstruct2);
            b = transform_poly(tmpcoeff2,tmpstruct2,polyb);
        end
    elseif any(strcmp(polya.type,{'monomial','Monomial'}))
        %                 | Domain A | Domain B   |
        %  ---------------------------------------|
        %   Type A (mono) |    ------------>      |
        %   -  -  -  -  - | -  -  -  | -  -  -  - |
        
        d = numel(a)-1;
        xa = polya.domain;
        xb = polyb.domain;
        z0 = (xb(1)*xa(2)-xb(2)*xa(1))/(xa(2)-xa(1));
        dz = (xb(2)-xb(1))/(xa(2)-xa(1));
        
        b = zeros(1,d+1);
        for m = 0:d
            for i = m:d
                b(end-m) = b(end-m) + a(end-i)/(dz)^i*nchoosek(i,i-m)*(-z0)^(i-m);
            end
        end
    else
        %                 | Domain A | Domain B   |
        %  ---------------------------------------|
        %   Type A        |    |     |      ^     |
        %   -  -  -  -  - | -  |  -  | -  - |-  - |
        %   Type B (mono) |    v----------->|     |
        %   -  -  -  -  - | -  -  -  | -  -  -  - |
        tmpstruct1 = struct('type','monomial','domain',polya.domain);
        tmpcoeff1 = transform_poly(a,polya,tmpstruct1);
        tmpstruct2 = struct('type','monomial','domain',polyb.domain);
        tmpcoeff2 = transform_poly(tmpcoeff1,tmpstruct1,tmpstruct2);
        b = transform_poly(tmpcoeff2,tmpstruct2,polyb);
    end
else
    %                 | Domain A |
    %  ---------------------------
    %   Type A        |    |     |
    %   -  -  -  -  - | -  |  -  |
    %   Type B        |    v     |
    %   -  -  -  -  - | -  -  -  |
    
    Pa = genpolyrecmatrix(numel(a)-1,polya.type);
    Pb = genpolyrecmatrix(numel(a)-1,polyb.type);
    b = (Pb\(Pa*a(:))).';
end

end

function P = genpolyrecmatrix(d,type)
% Generate the recursive transformation matrix from the given type to the
% monomial basis, given the degree d.

switch type
    case {'Monomial','monomial'}
        P = eye(d+1);
    case {'Legendre','legendre'}
        if d==0
            P = 1;
        else
            P = diag([1 1 zeros(1,d-1)]);
            for i = 3:d+1
                n = i-2;
                P(:,i) = [0;(2*n+1)/(n+1)*P(1:end-1,i-1)]-n/(n+1)*P(:,i-2);
            end
            P = P(end:-1:1,end:-1:1);
        end
    case {'Chebyshev','chebyshev'}
        if d==0
            P = 1;
        else
            P = diag([1 1 zeros(1,d-1)]);
            for i = 3:d+1
                P(:,i) = [0;2*P(1:end-1,i-1)]-P(:,i-2);
            end
            P = P(end:-1:1,end:-1:1);
        end
    case {'Chebyshev2','chebyshev2'}
        if d==0
            P = 1;
        else
            P = diag([1 2 zeros(1,d-1)]);
            for i = 3:d+1
                P(:,i) = [0;2*P(1:end-1,i-1)]-P(:,i-2);
            end
            P = P(end:-1:1,end:-1:1);
        end 
    otherwise
        error('transform_poly:genpolyrecmatrix:type',...
            'Unsupported type of polynomial basis');
end

end