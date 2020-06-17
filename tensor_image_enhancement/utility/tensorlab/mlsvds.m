function [U,S,sv] = mlsvds(T,size_core,largeScale,usesvds)
%MLSVDS Multilinear singular value decomposition for sparse tensors.
%   [U,S] = mlsvds(T, size_core) computes the factor matrices U{1}, ..., U{N}
%   and all-unitary and ordered core tensor S of size size_core belonging to a
%   truncated MLSVD of the sparse N-th order tensor T. The mode-n multilinear
%   singular values are stored in sv{n}. The vector sv{n} contains the Frobenius
%   norms of the mode-n slices of the core tensor S, e.g. S(:,:,k) are the
%   mode-3 slices of S.
%
%   [U,S,sv] = mlsvds(T, size_core) also computes the multilinear singular
%   values. The mode-n multilinear singular values are stored in sv{n}. The
%   vector sv{n} contains the Frobenius norms of the mode-n slices of the core
%   tensor S, e.g. S(:,:,k) are the mode-3 slices of S.
%
%   [U,S] = mlsvds(T, size_core, largeScale) can be used to disable (false) or
%   enable (true) the large scale implementation. If not set, the method
%   selects the largeScale option automatically depending on the size of the
%   tensor. 
% 
%   [U,S] = mlsvds(T, size_core, largeScale, usesvds) can be used disable
%   (false) or enable (true) the use of svds to compute the factor matrices
%   U{n}. Using svds can result in more accurate results, but is often (much)
%   slower. Default is false.
%
%   See also mlsvd, mlsvd_rsi, lmlra, tmprod, svd.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. De Lathauwer, B. De Moor, J. Vandewalle, "A Multilinear Singular
%       Value Decomposition", SIAM J. Matrix Anal. Appl., Vol. 21, No. 4,
%       April 2000, pp. 1253-1278.

    if isnumeric(T)
        T = fmt(T);
    end
    if ~isstruct(T) || ~isfield(T, 'sparse') || ~T.sparse
        error('mlsvds:T', ['T should be in the sparse tensor format (see ' ...
                           'fmt)']);
    else 
        T = fmt(T);
    end 
    N = length(T.size);
    if length(size_core) ~= N,
        error('mlsvds:size_core', ['length(size_core) should be ' ...
                            'length(T.size)']);
    end     
    if ~exist('largeScale', 'var') || ~isscalar(largeScale) || ...
            ~islogical(largeScale) 
        largeScale = length(T.val)*prod(size_core) > 1e9/8;
    end
    if ~exist('usesvds', 'var')
        usesvds = any(size_core == T.size);
    end

    %% Compute multilinear singular values and factor matrices
    U = cell(1,N);
    sv = cell(1,N);
    if usesvds 
        for n = 1:N 
            Tn = tens2mat(T, n);
            [U{n}, sv{n}, ~] = svds(Tn, size_core(n));
            sv{n} = diag(sv{n});
        end
    else 
        for n = 1:N 
            Tn = tens2mat(T, n);
            Tn = Tn*Tn';
            [U{n}, sv{n}] = eigs(Tn, size_core(n));
            sv{n} = sqrt(abs(diag(sv{n})));
            % eigs does not always sort the eigenvalues
            [sv{n}, i] = sort(sv{n}, 'descend');
            U{n} = U{n}(:,i);
        end
    end 

    %% Compute S
    % S is computed by multiplying T in each mode by U{n};
    if ~largeScale
        S = T.val;
        subs = T.sub;
        for n = N:-1:2
            S = bsxfun(@times, reshape(S,length(subs{1}), 1, prod(size_core(n+1:N))), ...
                       conj(U{n}(subs{n},:)));
            S = reshape(S, [], prod(size_core(n:N)));
            if n > 2
                ind = sub2ind(T.size(1:n-1), subs{1:n-1});
            else 
                ind = subs{1};
            end
            [~,ia,ic] = unique(ind);
            for k = 1:size(S,2)
                S(1:length(ia),k) = accumarray(ic, S(:,k));
            end
            
            S = S(1:length(ia),:);
            subs = cellfun(@(s) s(ia), subs(1:n-1), 'UniformOutput', false);
        end
        % n = 1
        S = bsxfun(@times, reshape(S,length(subs{1}),1, prod(size_core(2:N))), ...
                   conj(U{1}(subs{1},:)));
        S = reshape(sum(S), size_core);
    else
        %% Compute S memory efficient
        coregrid = arrayfun(@(n) 1:n, size_core, 'UniformOutput', false);
        [coregrid{:}] = ndgrid(coregrid{:});
        coregrid = cellfun(@(c) c(:), coregrid, 'UniformOutput', false);
        coregrid = cat(2, coregrid{:});
        caches = cell(1,N);
        S = zeros(size_core);

        % precompute the contraction indices
        subs = T.sub;
        ic = cell(1,N);
        for n = N:-1:2
            if n > 2
                ind = sub2ind(T.size(1:n-1), subs{1:n-1});
            else 
                ind = subs{1};
            end
            [~,ia,ic{n}] = unique(ind);
            for k = 1:n-1
                subs{k} = subs{k}(ia);
            end
        end

        prevgrid = -ones(1, N);
        caches{N} = T.val(:);
        for i = 1:prod(size_core)
            % update caches
            new = find(coregrid(i,:) ~= prevgrid, 1, 'last');
            for k = new-1:-1:1
                caches{k} = caches{k+1}.*conj(U{k+1}(subs{k+1},coregrid(i,k+1)));
                caches{k} = accumarray(ic{k+1}, caches{k});
            end
            S(i) = sum(caches{1}.*conj(U{1}(subs{1},coregrid(i,1))));
            
            prevgrid = coregrid(i, :);
        end
    end
end