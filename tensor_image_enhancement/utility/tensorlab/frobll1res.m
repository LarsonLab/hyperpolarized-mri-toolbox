function f = frobll1res(T, U, varargin)
%FROBLL1RES Frobenius norm of the residual of a LL1 decomposition.
%
%   FROBLL1RES(T,U,options), FROBLL1RES(T,U,L,options), FROBLL1RES(T,U, 'key',
%   value) or FROBLL1RES(T,U,L,'key',value) can be used to set the following
%   options:
%
%   - ExpandLimit           Expands structured tensors to full tensors if the
%                           total number of entries prod(getsize(T)) is
%                           smaller than the value of this limit. Default:
%                           1e6.
%   - Type                  The type of the tensor. If not given, the type is
%                           determined automatically. 
%
%   All other options are passed on to LL1RES if needed.
%
%   See also ll1res, frob, getstructure
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/30   NV      Initial version

    if nargin >= 3 && ~ischar(varargin{1})
        L = varargin(1);
        varargin = varargin(2:end);
    else 
        L = {};
    end
    
    p = inputParser;
    p.addOptional('ExpandLimit', 1e6);
    p.addOptional('Type', getstructure(T));
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    passopt = p.Unmatched;
    
    type = options.Type;
    if any(strcmpi(type, {'full', 'incomplete', 'sparse'}))
        f = frob(ll1res(T, U, L{:}, passopt));
    elseif prod(getsize(T)) <= options.ExpandLimit
        if ~isfield(passopt, 'Format'), passopt.Format = false; end
        f = frob(ll1res(ful(T), U, L{:}, passopt));
    else 
        if all(cellfun(@isnumeric, U))
            if isempty(L)
                error('frobll1res:L', ['U is given in the CPD format but no L ' ...
                                    'is given'])
            end
            if iscell(L), L = L{1}; end
            if length(U) ~= 3
                error('frobll1res:U', ...
                      'In the CPD format, U should have three factor matrices.');
            end 
            if any(cellfun('size', U(:).', 2) ~= [sum(L), sum(L), length(L)])
                error('frobll1res:U', ['cellfun(''size'', U, 2) should be [sum(L), ' ...
                                    'sum(L), length(L)])']);
            end
            % CPD format
            expvec = arrayfun(@(n) ones(1, L(n))*n, 1:length(L), 'UniformOutput',0); 
            expvec = cat(2, expvec{:});
            U{3} = U{3}(:,expvec);
            f = frob(T,'squared') - 2*real(inprod(T,U)) + frob(U,'squared');
        else 
            % BTD format
            f = frob(T,'squared') - 2*real(inprod(T,U)) + frob(U,'squared');
        end
        f = sqrt(abs(f));
    end
end
