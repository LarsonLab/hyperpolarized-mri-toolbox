function [U,S,sv] = mlsvd_rsi(T,size_core,varargin)
%MLSVD_RSI Sequentially truncated MLSVD using randomized subspace iteration.
%   [U,S] = MLSVD_RSI(T,size_core) computes the factor matrices U{1}, ..., U{N}
%   and all-unitary and ordered core tensor S of size size_core belonging to a
%   truncated MLSVD of the full or sparse tensor T, using sequential truncation
%   and using a randomized SVD algorithm based on randomized subspace
%   iteration.
%
%   [U,S,sv] = MLSVD_RSI(T,size_core) also computes the multilinear singular
%   values sv{n} for each mode n. The vector sv{n} contains the Frobenius norms
%   of the mode-n slices of the core tensor S, e.g. S(:,:,k) are the mode-3
%   slices of S.
%
%   [U,S] = MLSVD_RSI(T,size_core,options) or [U,S] = MLSVD_RSI(T, size_core,
%   key, value,...) can be used to set the following options:
%   - Compress = true   Remove the parts of the factor matrices U{n} and core
%                       tensor S corresponding due to the oversampling. 
%   - p = 5             The oversampling parameter. Either a scalar or a
%                       vector of length N = getorder(T).
%   - q = 2             Number of subspace iterations to be performed. Set to
%                       1 to reduce computation time (and accuracy!).

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Implementation note: the MLSVD [1] is computed using sequential trunction
%   [2] and randomized SVDs [3, alg 5.1] with randomized subspace iteration [3,
%   alg 4.4] for improved efficiency. 
%
%   References: 
%   [1] L. De Lathauwer, B. De Moor, J. Vandewalle, "A Multilinear Singular
%       Value Decomposition", SIAM J. Matrix Anal. Appl., Vol. 21, No. 4,
%       April 2000, pp. 1253-1278.
%   [2] N. Vannieuwenhoven, R. Vandebril, K. Meerbergen, "A new truncation
%       strategy for the higher-order singular value decomposition", SIAM
%       J. Sci. Comput., Vol. 34, No. 2, April 2012, pp. A1027-A1052.
%   [3] N. Halko, P.G. Martinsson, J.A. Tropp, "Finding Structure with
%       Randomness: Probabilistic Algorithms for Constructing Approximate Matrix
%       Decompositions", SIAM Rev, Vol.53, No.2, 2011, pp. 217-288.
    
    p = inputParser;
    p.addOptional('p', 5);
    p.addOptional('q', 2);
    p.addOptional('Compress', true);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    size_tens = getsize(T);
    N = length(size_tens);    
    
    if length(options.p) == 1
        p = ones(1,N)*options.p;
    elseif length(options.p) == N
        p = options.p(:).';
    else 
        error('mlsvd_rsi:invalidP', ...
              'p must be either a scalar or a vector of length getorder(T)');
    end
    q = options.q;
    
    % input checks
    if length(size_core) > N
        error('mlsvd_rsi:size_core', 'length(size_core) should be getorder(T)');
    end

    % check tensor type
    if ~isnumeric(T) && ~(isstruct(T) && isfield(T, 'sparse') && T.sparse)
        error('mlsvd_rsi:invalidType', ...
              'MLSVD_RSI can only be used for full or sparse tensors.');
    end
    Y = T;
    
    % outputs 
    U = cell(1, N);
    sv = cell(1,N);
    
    % compute compression
    for n = 1:N
        % randomized compression
        sz = [prod(size_tens([1:n-1 n+1:N])) size_core(n)+p(n)];
        Y = tens2mat(Y,n);
        tmp = Y*randn(sz);
        [Q,~] = qr(tmp,0);
        % subspace iteration
        for j = 1:q
            tmp = Y'*Q;
            [Q,~] = qr(tmp, 0);
            tmp = Y*Q;
            [Q,~] = qr(tmp, 0);
        end
        % SVD stage
        tmp = Q'*Y;
        [U{n},S,V] = svd(tmp, 'econ');
        U{n} = Q*U{n};
        sv{n} = diag(S);
        % prepare for next loop
        size_tens(n) = size(U{n},2);
        Y = S*V';
        Y = mat2tens(Y,size_tens,n);
        % compress if needed
        if options.Compress;
            U{n} = U{n}(:, 1:min(size_tens(n), size_core(n)));
        end
    end
    % compute core
    if ~options.Compress, size_core = size_core + p; end
    ind = arrayfun(@(n) 1:n, size_core(1:N), 'UniformOutput', false);
    S = Y(ind{:});
end
