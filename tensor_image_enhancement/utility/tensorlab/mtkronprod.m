function M = mtkronprod(T,U,n,transpose)
%MTKRONPROD Matricized tensor Kronecker product.
%   mtkronprod(T,U,n) computes the product
%
%      tens2mat(T,n)*conj(kron(U([end:-1:n+1 n-1:-1:1])))
%
%   without permuting the tensor T, or by exploiting the efficient
%   representation of structured tensors. Note that for improved performance, it
%   is advisable for the two largest dimensions of T to be the first and last
%   modes of T (in the case of full tensors).
%
%   mtkronprod(T,U,n,'T') and mtkronprod(T,U,n,'H') transpose or complex
%   conjugate transpose the matrices U{n} before computing the matricized
%   tensor Kronecker product.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2016/01/15   OD      Support structured tensors: Hankel
% - 2015/12/19   NV      Support structured tensors: CPD, BTD, LMLRA, TT

if nargin < 4, transpose = 0; end
if ischar(transpose)
    if strcmp(transpose,'T'), transpose = 1;
    elseif strcmp(transpose,'H'), transpose = 1i; end
end
U = U(:).';
if transpose == 0, size_tens = cellfun('size',U,1);
else size_tens = cellfun('size',U,2); end
data = T;
isstructured = false;
switch getstructure(T)
    case {'cpd','btd','lmlra','tt','hankel','loewner'}
        isstructured = true;
    case {'incomplete','sparse'}
        T = data.matrix;
        alloc = @sparse;
    case 'full'
        alloc = @zeros;
    otherwise
        error('mtkrprod:notImplemented', 'This is not yet implemented');
end

N = length(size_tens);
if transpose == 0, size_core = cellfun('size',U,2);
else size_core = cellfun('size',U,1); end
if n > 0, size_core(n) = size_tens(n); end
R = prod(size_core([1:n-1 n+1:N]));

if isstructured
    type = getstructure(T);
    idx = [1:n-1 n+1:length(T)];
    switch type
        case 'cpd'
            switch transpose
                case 0
                    W = cellfun(@(u,v) v'*u, T(idx), U(idx), 'UniformOutput', ...
                        false);
                case 1
                    W = cellfun(@(u,v) conj(v)*u, T(idx), U(idx), 'UniformOutput', ...
                        false);
                case 1i
                    W = cellfun(@(u,v) v*u, T(idx), U(idx), 'UniformOutput', ...
                        false);
            end
            W = kr(W(end:-1:1));
            if n > 0, M = T{n}*W.';
            else M = sum(W,2).'; end
        case 'lmlra'
            W = cell(1,length(T{1}));
            for m = 1:length(T{1})
                if m == n, W{n} = T{1}{n};
                else
                    switch transpose
                        case 0,  W{m} = U{m}'*T{1}{m};
                        case 1,  W{m} = conj(U{m})*T{1}{m};
                        case 1i, W{m} = U{m}*T{1}{m};
                    end
                end
            end
            M = mtkronprod(T{2},W,n,'H');
            if n > 0, M = T{1}{n}*M; end
        case 'btd'
            M = 0;
            for r = 1:length(T)
                M = M + mtkronprod([{T{r}(1:end-1)}, T{r}{end}], U, n, transpose);
            end
        case 'tt'
            ttr = [1 cellfun('size', T(2:end), 1) ones(1, length(U)-length(T)+1)];
            if ~isreal(transpose),
                U = cellfun(@conj, U, 'UniformOutput', false);
            end
            if transpose ~= 0
                U = cellfun(@(u) u.', U, 'UniformOutput', false);
            end
            
            tmp = cell(1, length(U));
            if n ~= 1, tmp{1} = U{1}'*T{1};
            else tmp{1} = T{1}; end
            for m = 2:length(T)
                if m == n, tmp{n} = T{n}; continue; end
                tmp{m} = reshape(permute(T{m}, [2 1 3]), size_tens(m), ttr(m)*ttr(m+1));
                tmp{m} = reshape(U{m}'*tmp{m},size_core(m),ttr(m),ttr(m+1));
                tmp{m} = permute(tmp{m}, [2 1 3]);
            end
            M = ttgen(tmp);
            if n > 0, M = tens2mat(M, n);
            else M = M(:).'; end
        case 'hankel'
            HN = size(T.val,1);
            U = cellfun(@conj,U,'UniformOutput',false); % needed to turn
            if ~isreal(transpose),
                U = cellfun(@conj, U, 'UniformOutput', false);
            end
            if transpose ~= 0
                U = cellfun(@(u) u.', U, 'UniformOutput', false);
            end
            if T.ispermuted
                others = T.order+(1:numel(T.subsize.other));
                if n<=T.order
                    % T is permuted and n is in 1:T.order
                    Tval = mat2tens(T.val,[size(T.val,1),T.subsize.other],1);
                    
                    % Using Khatri-Rao product for other dimensions
                    Tkr = tmprod(Tval,U(others),2:ndims(Tval),'T');
                    Tmat = tens2mat(Tkr,1);
                    
                    % Using FFT for tensorized dimensions
                    F = fft(Tmat,HN,1);
                    UU = kr(cellfun(@(x) fft(x(end:-1:1,:),HN,1).',U([T.order:-1:n+1 n-1:-1:1]),'UniformOutput',false));
                    TT = kr(F.',UU);
                    F = ifft(TT.',[],1);
                    if n~=0, M = F(end-size(U{n},1)+1:end,:);
                    else M = F(end,:); end
                else
                    % T is permuted and n is larger than T.order
                    othersn = [T.order+1:n-1 n+1:numel(T.size)];
                    
                    % Using FFT for tensorized dimensions
                    F = fft(U{1},HN,1).';
                    for i = 2:T.order
                        F = kr(fft(U{i},HN,1).',F);
                    end
                    F = F.';
                    tmp = ifft(F);
                    tmp = T.val.'*tmp;
                    if ~isempty(othersn)
                        tmpT = mat2tens(tmp.',[size(tmp,2) T.subsize.other],1);
                        % Using Khatri-Rao product for other dimensions
                        M = tmprod(tmpT,U(othersn),[2:n-T.order n-T.order+2:ndims(tmpT)],'T');
                        M = tens2mat(M,n-T.order+1);
                    else
                        if n==0, M = reshape(tmp,1,[]);
                        else M = tmp; end
                    end
                end
            else
                others = [1:T.dim-1 (T.dim+T.order):numel(T.size)];
                if n>=T.dim && n<T.dim+T.order
                    if ~isempty(others)
                        Tval = mat2tens(T.val,[size(T.val,1),T.subsize.other],1);
                        Tkr = tmprod(Tval,U(others),2:ndims(Tval),'T');
                        Tmat = tens2mat(Tkr,1);
                    else
                        Tmat = T.val;
                    end
                    vec = [T.dim:n-1 n+1:T.dim+T.order-1];
                    
                    % Using Khatri-Rao product for other dimensions
                    UU = kr(cellfun(@(x) fft(x(end:-1:1,:),HN,1).',U(vec(end:-1:1)),'UniformOutput',false));
                    
                    % Using FFT for tensorized dimensions
                    F = fft(Tmat,HN,1).';
                    TT = kr(F,UU).';
                    F = ifft(TT,[],1);
                    M = F(end-size(U{n},1)+1:end,:);
                    
                    M = mat2tens(M,[size(M,1) cellfun('size',U([T.dim:n-1 n+1:T.dim+T.order-1 others]),2)],1);
                    vec = [1 T.order+1:T.dim+T.order-1 2:T.order T.dim+T.order:ndims(M)];
                    M = permute(M,vec);
                    M = tens2mat(M,1);
                else
                    if n<T.dim
                        % T is not permuted and n is smaller than the Hankel dimensions
                        othersn = [1:n-1 n+1:T.dim-1 (T.dim+T.order):numel(T.size)];
                        vec = [1 T.order+2:T.order+T.dim-1 2:T.order+1 T.order+T.dim:numel(T.size)];
                        tmprodvec = [2:n n+2:numel(T.subsize.other)+1];
                        tens2matidx = n+1;
                    else % n>=T.dim+T.order
                        % T is not permuted and n is larger than the Hankel dimensions
                        othersn = [1:T.dim-1 (T.dim+T.order):n-1 n+1:numel(T.size)];
                        vec = [1 T.order+2:T.order+T.dim 2:T.order+1 T.order+T.dim+1:numel(T.size)];
                        tmprodvec = [2:n-T.order n-T.order+2:numel(T.subsize.other)+1];
                        tens2matidx = n-T.order+1;
                    end
                    
                    % Using FFT for tensorized dimensions
                    F = fft(U{T.dim},HN,1).';
                    for i = 1:T.order-1
                        F = kr(fft(U{T.dim+i},HN,1).',F);
                    end
                    tmp = ifft(F.');
                    tmp = tmp.'*T.val;
                    if ~isempty(othersn)
                        M = mat2tens(tmp,[size(tmp,1) T.subsize.other],1);
                        % Using Khatri-Rao product for other dimensions
                        M = tmprod(M,U(othersn),tmprodvec,'T');
                        if n~=0
                            M = tens2mat(M,tens2matidx);
                            M = mat2tens(M,[size(M,1) cellfun('size',U([T.dim:T.dim+T.order-1 othersn]),2)],1);
                            M = permute(M,vec);
                            M = tens2mat(M,1);
                        else
                            M = reshape(M,cellfun('size',U([T.dim:T.dim+T.order-1 othersn]),2));
                            M = ipermute(M,[T.dim:T.dim+T.order-1 1:T.dim-1 T.dim+T.order:numel(T.size)]);
                            M = reshape(M,1,[]);
                        end
                    else
                        if n==0, M = reshape(tmp,1,[]);
                        else M = tmp; end
                    end
                end
            end
        case 'loewner'
            error('mtkronprod:loewner','Loewner structure is not supported yet; use the ful version instead.');
        otherwise
            error('Not implemented yet!');
    end
    return
end

% Determine if large-scale version of the algorithm should be executed.
ratio = size_core./size_tens;
perm = zeros(1,N);
l = 1; r = N;
for i = 1:N
    if r == n || (l < n && ratio(l) < ratio(r)), perm(i) = l; l = l+1;
    else perm(i) = r; r = r-1;
    end
end
if n == 0, mem = 8*prod(size_tens);
else mem = 8*prod(size_tens([1:perm(1)-1 perm(1)+1:N]))*size_core(perm(1));
end
largescale = mem > 2e9 || (isstruct(data) && isempty(data.matrix));
if largescale && ~isstruct(data)
    error('mtkronprod:mem',['Intermediate result is too large. Try ' ...
        'to format data with fmt first so that the large-scale ' ...
        'version of the algorithm can be applied.']);
end

if largescale
    
    % The large scale implementation fires if the required memory would
    % otherwise be larger than 2GB and assumes the dataset has been
    % formatted by fmt.
    idx = [1:n-1 n+1:N];
    jdx = cell(1,N);
    switch transpose
        case 0,  U = cellfun(@(u)conj(u),U,'UniformOutput',false);
        case 1,  U = cellfun(@(u)u',U,'UniformOutput',false);
        case 1i, U = cellfun(@(u)u.',U,'UniformOutput',false);
    end
    if n == 0, M = zeros(1,R);
    else M = zeros(size_tens(n),R); end
    for r = 1:R
        s = data.val;
        [jdx{:}] = ind2sub(size_core([1:n-1 n+1:N]),r);
        for j = 1:length(idx)
            s = s.*U{idx(j)}(data.sub{idx(j)},jdx{j});
        end
        if n == 0, M(r) = sum(s);
        else M(:,r) = accumarray(double(data.sub{n}),s,[size_tens(n) 1]);
        end
    end
    
else
    
    % Apply structure-exploiting matricized tensor Kronecker product.
    M = T;
    cpl = cumprod([1 size_core(perm(perm < n))]); l = 1;
    cpr = cumprod([1 size_core(perm(perm > n))]); r = 1;
    for i = 1:length(perm)
        
        mode = perm(i);
        if mode < n
            
            tmp = reshape(M,cpl(l)*size_tens(mode),[]);
            M = alloc(cpl(l),size(tmp,2)*size_core(mode));
            for j = 1:cpl(l)
                idx = j:cpl(l):size(tmp,1);
                switch transpose
                    case 0,  tmp2 = U{mode}'*tmp(idx,:);
                    case 1,  tmp2 = conj(U{mode})*tmp(idx,:);
                    case 1i, tmp2 = U{mode}*tmp(idx,:);
                end
                M(j,:) = tmp2(:);
            end
            l = l+1;
            
        elseif mode > n
            
            tmp = reshape(M,[],cpr(r)*size_tens(mode));
            M = alloc(size(tmp,1)*size_core(mode),cpr(r));
            for j = 1:cpr(r)
                idx = (1:size_tens(mode))+(j-1)*size_tens(mode);
                switch transpose
                    case 0,  tmp2 = tmp(:,idx)*conj(U{mode});
                    case 1,  tmp2 = tmp(:,idx)*U{mode}';
                    case 1i, tmp2 = tmp(:,idx)*U{mode}.';
                end
                M(:,j) = tmp2(:);
            end
            r = r+1;
            
        end
        
    end
    
    % Permute and reshape output.
    if n > 1
        if issparse(M)
            idx = find(M);
            [i,j,k] = ind2sub([cpl(end) size_tens(n) cpr(end)],idx);
            ik = sub2ind([cpl(end) cpr(end)],i,k);
            M = sparse(j,ik,full(M(idx)),size_tens(n),cpl(end)*cpr(end));
        else
            M = reshape(M,[cpl(end) size_tens(n) cpr(end)]);
            M = reshape(permute(M,[2 1 3]),size_tens(n),[]);
        end
    else
        if n == 0, M = reshape(M,1,[]);
        else M = reshape(M,size_tens(n),[]); end
    end
    
end

end
