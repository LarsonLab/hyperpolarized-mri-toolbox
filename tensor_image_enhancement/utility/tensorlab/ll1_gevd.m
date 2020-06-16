function [U, output] = ll1_gevd(T, L, varargin)
%LL1_GEVD LL1 by generalized eigenvalue decomposition
%   U = LL1_GEVD(T,L) computes the (L,L,1)-decomposition of a full or sparse
%   tensor T. L is an array of length R. The result U is a cell of length R,
%   each containing a (L,L,1)-term. A (L,L,1)-term stores the factor matrices
%   A_r and B_r, and the factor vector C_r. A fourth identity matrix is added to
%   fit in the general BTD format.
%
%   The matrices A_r and B_r are grouped using an heuristic based on k-means,
%   which is nondeterministic. If the heuristic fails to group the matrices,
%   options.Backup is used as backup strategy.
%
%   [U, output] = LL1_GEVD(...) also returns an output struct containing the
%   following fields:
%
%      output.Backup          - True if options.Backup is used.
%
%   LL1_GEVD(T, L, options) or LL1_GEVD(T,L,'key',value,...) may be set to the
%   following options:
%
%      options.gevd =         - The GEVD solver used. 
%      [@eig|{@gevd_bal}]
%      options.maxRetries =   - The number of retries for the clustering
%      100
%      options.Backup =       - Back up algorithm to use when no LL1_GEVD
%      @ll1_rnd                 failed.
%      options.isReal =       - Compute a complex solution if isReal=false
%      [true|false]             even if T is real. If not set, isReal =
%                               isreal(T)
%      options.Slices =       - Method for computing the slices used in the
%      [{'mlsvd'}|'random']     GEVD. 'random' is recommended if the first two 
%                               frontal slices of the mlsvd of the (permuted)
%                               tensor T are have rank < sum(L).
%      options.OutputFormat = - Format of the result U.
%      [{'btd'}|'cpd']
    
%   Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%              Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%  
%   References:
%   [1] Lieven De Lathauwer, "Decompositions of a Higher-Order Tensor in Block
%       Terms --- Part II: Definitions and Uniqueness", SIAM J. Matrix Anal.
%       Appl. Vol. 30, No. 3, pp 1033-1066
    
% Version History:
% - 2016/03/20   NV      Random slices
% - 2016/01/04   NV      Complex tensors & better clustering
% - 2014/11/24   NV      Initial version
    
    p = inputParser;
    p.addOptional('gevd', @gevd_bal);
    p.addOptional('maxRetries', 100);
    p.addOptional('Backup', @ll1_rnd);
    p.addOptional('OutputFormat', 'btd');
    p.addOptional('Slices', 'mlsvd');
    p.addOptional('isReal', isreal(ful(T,1)));
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    R = length(L);
    sL = sum(L);
    
    if ~any(strcmpi(options.OutputFormat, {'cpd', 'btd'}))
        error('ll1_gevd:OutputFormat', 'Unknown output format %s.', ...
              options.OutputFormat);
    end 
    
    % Check some prerequisites
    if ~isnumeric(T) && ~(isstruct(T) && isfield(T, 'sparse') && T.sparse)
        error('This method currently only supports full and sparse tensors.');
    end
    if isnumeric(T), size_tens = size(T);
    else size_tens = T.size; end
    N = length(size_tens);

    if any(size_tens(1:2) < [sL, sL])
        error('ll1_gevd:T', 'size(T,n) should be >= sum(L) for n=1,2.');
    end
    if R == 1
        if isstruct(T)
            [U,S] = mlsvds(T, [sL,sL,1]);
        else
            [U,S] = mlsvd(T, [sL,sL,1]);
        end
        U{1} = U{1}*sqrt(S);
        U{2} = U{2}*sqrt(S);
        U = {[U, eye(sL)]};
        output.Backup = false;
        return;
    end
    if N ~= 3
        error('ll1_gevd:T', 'ndims(T) should be 3');
    end
    if size_tens(3) < 2
        error('ll1_gevd:T', 'size(T,3) should be >= 2.');
    end 
    
    % compress
    if strcmpi(options.Slices, 'random')
        size_core = [sL, sL, size_tens(3)];
    else 
        size_core = [sL, sL, min(max(R,2),size_tens(3))];
    end
        
    if isstruct(T), [Um, Sm, ~] = mlsvds(T, size_core);
    else [Um, Sm, ~] = mlsvd(T, size_core); end

    if strcmpi(options.Slices, 'random')
        Sm = tmprod(Sm, Um{3}.'*randn(size(Um{3},1), max(2,R)), 3, 'T');
    end

    % Compute GEVD to recover A, use D to make everything real if necessary
    [V,D] = options.gevd(Sm(:,:,1).', Sm(:,:,2).');
    D = diag(D);
    [~, i] = sort(sign(real(D)).*abs(D), 'descend');
    D = D(i);
    V = V(:, i);    
    if options.isReal && ~isreal(V)
        c = find(imag(D) ~= 0);
        for k = 1:2:length(c)
            V(:,c(k)) = real(V(:,c(k)));
            V(:,c(k+1)) = imag(V(:,c(k+1)));
            D(c(k+1)) = conj(D(c(k+1)));
        end
    end
    
    % Construct A and B
    A = pinv(V).';
    B = Sm(:,:,1).'*V;
    
    % Cluster vectors 
    if sum(L) ~= length(L)
        Ctmp = (pinv(kr(B,A))*tens2mat(Sm,[1 2], 3)).';
        Ctmp = bsxfun(@rdivide, Ctmp, sqrt(sum(abs(Ctmp).^2)));
        Ctmp(isnan(Ctmp)) = 0;
        [u,~,~] = svd(Ctmp,'econ');
        tmp = abs(u(:,1:R)'*Ctmp);
        for k = 1:options.maxRetries
            [l,~] = kmeans(tmp, R);
            L2 = histc(l, 1:R);
            if all(sort(L2)==sort(L)); break; end
        end 
    else 
        L2 = L;
        l = 1:sum(L);
    end

    output.Backup = false;
    if all(sort(L2) == sort(L))
        % compute C
        P = sparse(1:sL,l,1,sL,R);
        if ~strcmpi(options.Slices, 'random')
            C = (pinv(kr(B,A)*P)*tens2mat(Sm,[1 2],3)).';
            C = Um{3}*C;
        end

        % expand
        A = Um{1}*A;
        B = Um{2}*B;
        if strcmpi(options.Slices, 'random')
            C = (pinv(kr(B,A)*P)*tens2mat(T,[1 2],3)).';
        end
        
        % Format the output
        U = cell(1, R);
        for r = 1:R
            U{r}{1} = A(:,l==r);
            U{r}{2} = B(:,l==r);
            U{r}{3} = C(:,r);
            U{r}{4} = eye(sum(l==r));
        end

        % Sort terms to match given L
        [~, i] = sort(L, 'descend');
        [~, j] = sort(L2, 'descend');
        U(i) = U(j);
        
        if strcmpi(options.OutputFormat, 'cpd')
            U = ll1convert(U);
        end
    else 
        % could not cluster vectors, use backup
        U = options.Backup(T, L, 'OutputFormat', options.OutputFormat);
        output.Backup = true;
    end
end 
