function E = ll1res(T,U,varargin)
%LL1RES Residual for a LL1 decomposition.
%   E = LL1RES(T,U) computes the residual tensor E as ll1gen(U)-T, in which U is
%   given in the BTD format. If T is an incomplete tensor, the residual is only
%   computed for known elements. If T is an efficient representation of a
%   structured tensor, T is first expanded using ful(T).
%
%   E = LL1RES(T,U,L) computes the residual tensor E as ll1gen(U,L)-T, in which U
%   is given in the CPD format.
%
%   LL1RES(T,U,options), LL1RES(T,U,L,options), LL1RES(T,U,'key',value) and
%   LL1RES(T,U,L,'key',value) can be used to set the following options:
%   
%     options.Format = [{true},false]   - If true, the tensor is formatted
%                                         using fmt(T) before computing the
%                                         residual.
%     options.Type                      - The type of the tensor. If not
%                                         given, the type is determined
%                                         automatically
%     options.TolLargeScale = 0.02      - If length(T.ind)/numel(T) <=
%                                         TolLargeScale, a large scale
%                                         variant of the algorithm is used.
%   See also cpdres, btdres, lmlrares
    
%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version History:
%   - 2016/01/05   NV      Initial version
    
    if nargin >= 3 && ~ischar(varargin{1}) && ~isstruct(varargin{1})
        L = varargin{1};
        varargin = varargin(2:end);
    else 
        L = [];
    end
    
    if all(cellfun(@isnumeric, U)) 
        if isempty(L)
            error('ll1res:L', ...
                  'If U is given in the CPD format, L should be given');
        end
        if length(U) ~= 3
            error('ll1res:U', ...
                  'In the CPD format, U should have three factor matrices.');
        end 
        if any(cellfun('size', U(:).', 2) ~= [sum(L), sum(L), length(L)])
            error('ll1res:U', ['cellfun(''size'', U, 2) should be [sum(L), ' ...
                               'sum(L), length(L)])']);
        end
        expvec = arrayfun(@(n) ones(1, L(n))*n, 1:length(L), 'UniformOutput', false); 
        expvec = cat(2, expvec{:});
        U{3} = U{3}(:,expvec);
        E = cpdres(T,U,varargin{:});
    else 
        E = btdres(T,U,varargin{:});
    end
end
