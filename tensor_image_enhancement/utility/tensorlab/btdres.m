function E = btdres(T,U,varargin)
%BTDRES Residual of a BTD.
%   E = BTDRES(T,U) computes the residual tensor E as btdgen(U)-T. If T is an
%   incomplete tensor, the residual is only computed for known elements. If T
%   is an efficient representation of a structured tensor, T is first expanded
%   using ful(T).
%
%   BTDRES(T, U, options) can be used to set the following options:
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
%
%   See also cpdres, lmlrares, ll1res, frobbtdres.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
p = inputParser;
p.addOptional('TolLargeScale', 0.02);
p.addOptional('Format', true);
p.addOptional('Type', getstructure(T));
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;

% Compute the residual.
if any(strcmpi(options.Type, {'full', 'incomplete', 'sparse'})) 
    if options.Format, T = fmt(T,true); end
else 
    T = ful(T);
end
if isstruct(T)
    if ~T.incomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        E = U{1}{1}*mtkronprod(U{1}{end},U{1}(1:end-1),1,'H');
        for r = 2:length(U)
            E = E+U{r}{1}*mtkronprod(U{r}{end},U{r}(1:end-1),1,'H');
        end
        if T.incomplete, E = E(T.ind)-T.val;
        elseif T.sparse, 
            E(T.ind) = E(T.ind)-T.val;
            E = reshape(E, T.size);
        else E = E-reshape(T,size(E));
        end
    else
        E = -T.val;
        for r = 1:length(U)
            S = U{r}{end};
            size_core = cellfun('size',U{r}(1:end-1),2);
            idx = cell(1,length(U{r})-1);
            for i = 1:numel(S)
                [idx{:}] = ind2sub(size_core,i);
                tmp = S(idx{:})*U{r}{1}(T.sub{1},idx{1});
                for n = 2:length(size_core)
                    tmp = tmp.*U{r}{n}(T.sub{n},idx{n});
                end
                E = E+tmp;
            end
        end
    end
    if T.incomplete
        E = struct('matrix', T.matrix, 'size', T.size, 'ind', T.ind, 'sub', {T.sub}, ...
                   'val', E, 'incomplete', true, 'sparse', false);
        if ~isempty(T.matrix)
            E.matrix = sparse(double(E.sub{1}), ...
                              double(1+idivide(E.ind-1,int64(size(T.matrix,1)))), ...
                              E.val,size(T.matrix,1),size(T.matrix,2));
        end
    end
else
    E = btdgen(U)-T;
end
