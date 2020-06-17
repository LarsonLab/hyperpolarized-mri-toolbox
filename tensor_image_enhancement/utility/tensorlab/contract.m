function T = contract(T, U, modes)
%CONTRACT Mode-n tensor vector contraction.
%   S = CONTRACT(T, U, modes) computes the contraction of a full or sparse
%   tensor T with the vectors U{1}, ..., U{N} along modes mode(1), ..., mode(N)
%   for N = length(T). modes should be a vector with distinct integers between 1
%   and getorder(T). The result S is a tensor of order getorder(T) - length(U).
%   
%   See also tmprod.
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@esat.kuleuven.be)
%
% Version History:
% - 2014/12/22   NV      Initial version

    try     
        type = getstructure(T);
    catch 
        error('contract:unknownType', 'Could not determine the type of T');
    end
    if ~any(strcmpi(type, {'full', 'sparse'}))
        error('contract:notImplemented', ['Contract only works on full and ' ...
                            'sparse tensors']);
    end
    sz = getsize(T);
    
    N = length(sz);
    [modes, ~, i] = unique(modes);
    U(i) = U;
    U = cellfun(@(u) u(:).', U, 'UniformOutput', false);
    if any(modes > N) || any(modes < 1)
        error('contract:modes', ['modes should be integer and 1 <= modes(n) <=' ...
                            ' ndims(T) for all n=1,...,length(modes)']);
    end
    if any(cellfun(@length, U) ~= sz(modes))
        error('contract:U', 'length(U{n}) should be size(T, modes(n))');
    end
        
    if length(modes) < N - 2 && strcmp(type, 'full')
        T = tmprod(T, U, modes, 0, 'saveperm');
        T = squeeze(T);
        return;
    end
    
    if strcmp(type, 'sparse') 
        S = T.val(:);
        subs = cellfun(@(s) double(s(:)), T.sub, 'UniformOutput', false);
        subs = cat(2, subs{:});
        remainingmodes = true(1, N);
        for n = 1:length(modes) 
            tmp = find(remainingmodes);
            remainingmodes(modes(n)) = false;
            S = S .* U{n}(subs(:,tmp==modes(n))).';  
            subs = subs(:,tmp~=modes(n));
            if sum(remainingmodes) >= 2
                csz = cumprod([1 sz(remainingmodes)]);             
                % -1 below not necessary
                ind = sum(bsxfun(@times, subs-1, csz(1:end-1)),2); %-1
            elseif sum(remainingmodes) == 1
                ind = subs;
            else 
                ind = ones(1, length(S));
            end
            [ia,ic] = order(ind);
            S = accumarray(ic, S);
            subs = subs(ia,:);
        end
        if sum(remainingmodes) == 0
            T = S;
        elseif sum(remainingmodes) == 1
            T = zeros(sz(remainingmodes), 1);
            T(subs) = S;
        elseif sum(remainingmodes) == 2
            remainingmodes = find(remainingmodes);
            i = subs(:,1);
            j = subs(:,2);
            m = sz(remainingmodes(1));
            n = sz(remainingmodes(2));
            T = sparse(i,j, S, m, n);
            T = full(T);
        else 
            T = struct;
            T.size = sz(remainingmodes);
            T.sub = num2cell(subs,1);
            T.val = S;
            T.incomplete = 0;
            T.sparse = 1;
            T = fmt(T);
        end 
    else 
        bndleft = find(modes == 1:length(modes), 1, 'last');
        bndright = find(modes == (N-length(modes)+1):N, 1, 'first');
        if isempty(bndleft), bndleft = 0; end
        if bndright < bndleft, bndright = bndleft + 1; end
        
        % Left-to-right contractions
        for n = 1:bndleft
            T = reshape(T, sz(n), prod(sz(n+1:end)));
            T = U{n}*T;        
        end
        sz(1:bndleft) = 1;
        
        % Right-to-left contractions
        for n = length(modes):-1:bndright
            T = reshape(T, prod(sz(1:modes(n)-1)), sz(modes(n)));
            T = T*U{n}.';        
        end
        sz(modes(bndright:end)) = [];
        sz(modes(1:bndleft)) = [];
        
        % Middle contractions from left to right
        if length(sz) >= 2
            T = reshape(T, sz);
            T = permute(T, [2:length(sz)-1, 1, length(sz)]);
            sz = size(T);

            for n = 1:length(sz)-2
                T = reshape(T, sz(n), prod(sz(n+1:end)));
                T = U{n+bndleft}*T;  
            end
        end
        if length(sz) > 1, T = reshape(T, sz(end-1:end)); end
    end
    
    function [ia,ic] = order(a)
        if issorted(a)
            a = a(:);
            indSortA = (1:length(a))';
        else 
            [a,indSortA] = sort(a(:));
        end 
        groupsSortA = [true; diff(a) ~= 0];
        ia = indSortA(groupsSortA);
        ic = cumsum(groupsSortA);
        ic(indSortA) = ic;
    end
    
end
