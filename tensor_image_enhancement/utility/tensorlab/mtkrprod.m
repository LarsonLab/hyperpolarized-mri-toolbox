function M = mtkrprod(T,U,n,conjugate)
%MTKRPROD Matricized tensor Khatri-Rao product.
%   mtkrprod(T,U,n) computes the product
%
%      tens2mat(T,n)*conj(kr(U([end:-1:n+1 n-1:-1:1])))
%
%   without permuting the tensor T, or by exploiting the efficient
%   representation of structured tensors.. Note that for improved performance,
%   it is advisable for the two largest dimensions of T to be the first and last
%   modes of T (in the case of full tensors). 
%
%   mtkrprod(T,U,n,false) does not apply the complex conjugate to the
%   factor matrices U{n} before computing the matricized tensor Khatri-Rao
%   product.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] N. Vannieuwenhoven, N. Vanbaelen, K. Meerbergen, R. Vandebril,
%       "The dense multiple-vector tensor-vector product: An initial
%       study," Technical Report TW635, Dept. of Computer Science,
%       KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2016/01/14   OD      Support structured tensors: Hankel, Loewner
% - 2015/12/19   NV      Support structured tensors: CPD, BTD, LMLRA, TT

if nargin < 4, conjugate = true; end
data = T;
size_tens = cellfun('size',U,1);
isstructured = false;


N = length(size_tens);
idx = [1:n-1 n+1:N];
R = size(U{idx(1)},2);

switch getstructure(T)
    case {'cpd','btd','lmlra','tt'}
        isstructured = true;
        nwas0 = n == 0;
        if nwas0, n = 1; idx = [1:n-1 n+1:N]; end
    case {'hankel','loewner'}
        isstructured = true;
        nwas0 = false;
    case {'incomplete','sparse'}
        nwas0 = n == 0;
        if nwas0, n = 1; idx = [1:n-1 n+1:N]; end
        T = data.matrix;
        alloc = @sparse;
    case 'full'
        alloc = @zeros;
    otherwise
        error('mtkrprod:notImplemented', 'This is not yet implemented');
end
if isstructured
    switch getstructure(T)
        case 'cpd'
            for m = length(T)+1:length(U)
                T{m} = ones(1,size(T{1},2));
            end
            if conjugate
                W = cellfun(@(u,v) u'*v, T(idx), U(idx), 'UniformOutput', false);
                W = prod(cat(3, W{:}),3);
                W = conj(W);
            else
                W = cellfun(@(u,v) u.'*v, T(idx),U(idx), 'UniformOutput', false);
                W = prod(cat(3, W{:}),3);
            end
            M = T{n}*W;
        case 'lmlra'
            W = cell(1, length(T{1}));
            for m = 1:length(T{1})
                if m == n
                    W{m} = T{1}{n};
                    continue
                end
                if conjugate, W{m} = T{1}{m}'*U{m};
                else W{m} = T{1}{m}.'*U{m}; end
            end
            W = [W U(length(T{1})+1:end)];
            M = mtkrprod(T{2}, W, n, conjugate);
            if n <= length(T{1})
                M = T{1}{n}*M;
            end
        case 'btd'
            M = 0;
            for r = 1:length(T)
                M = M + mtkrprod([{T{r}(1:end-1)}, T{r}{end}], U, n, conjugate);
            end
        case 'tt'
            ttr = [1 cellfun('size', T(2:end), 1) ones(1, length(U)-length(T)+1)];
            size_tens = cellfun('size', U, 1);
            
            if conjugate, U = cellfun(@conj, U, 'UniformOutput', false); end
            
            tmp = cell(1, length(U));
            if n ~= 1, tmp{1} = T{1}.'*U{1}; end
            for m = 2:length(T)
                if m == n, continue; end
                tmp{m} = reshape(permute(T{m}, [1 3 2]), ttr(m)*ttr(m+1), size_tens(m));
                tmp{m} = reshape(tmp{m}*U{m},ttr(m),ttr(m+1),R);
            end
            for m = length(T)+1:length(U)
                tmp{m} = reshape(U{m}, 1, 1, R);
                T{m} = 1;
            end
            M = nan(size(U{n}));
            for r = 1:R
                left = 1;
                if n ~= 1, left = tmp{1}(:,r).'; end
                for m = 2:n-1
                    left = left * tmp{m}(:,:,r);
                end
                right = 1;
                for m = n+1:length(U)
                    right = right * tmp{m}(:,:,r);
                end
                left = left*reshape(T{n},ttr(n),size_tens(n)*ttr(n+1));
                M(:,r) = reshape(left, size_tens(n), ttr(n+1))*right;
            end
        case 'hankel'
            HN = size(T.val,1);
            if conjugate, U = cellfun(@(x) conj(x),U,'UniformOutput',false); end
            
            if T.ispermuted
                others = T.order+(1:numel(T.subsize.other));
                if n<=T.order
                    % T is permuted and n is in 1:T.order
                    F = fft(T.val*kr(U{others(end:-1:1)}),HN,1);
                    for i = [1:n-1 n+1:T.order]
                        F = F.*fft(U{i}(end:-1:1,:),HN,1);
                    end
                    F = ifft(F,[],1);
                    if n~=0, M = F(end-size(U{n},1)+1:end,:);
                    else M = F(end,:); end
                else
                    % T is permuted and n is larger than T.order
                    othersn = [T.order+1:n-1 n+1:numel(T.size)];
                    F = ones(HN,R);
                    for i = 1:T.order
                        F = F.*fft(U{i},HN,1);
                    end
                    tmp = ifft(F);
                    tmp = T.val.'*tmp;
                    if ~isempty(othersn)
                        % If there are other modes besides the tensorized modes
                        M = zeros(size(U{n}));
                        for r = 1:R
                            newT = reshape(tmp(:,r),T.subsize.other);
                            ur = cellfun(@(x) x(:,r),U(othersn),'UniformOutput',false);
                            M(:,r) = contract(newT,ur,[1:n-T.order-1 n-T.order+1:ndims(newT)]);
                        end
                    else M = tmp;
                    end
                end
            else
                others = [1:T.dim-1 (T.dim+T.order):numel(T.size)];
                if n>=T.dim && n<T.dim+T.order
                    % T is not permuted and n is inside the Hankel modes
                    if ~isempty(others)
                        % If there are other modes besides the tensorized modes
                        tmp = T.val*kr(U{others(end:-1:1)});
                    else
                        % tmp = T.val*ones(size(T.val,2),1);
                        tmp = repmat(T.val,1,R);
                    end
                    F = fft(tmp,HN,1); % TODO: mtkrprod
                    for i = [T.dim:n-1 n+1:(T.dim+T.order-1)]
                        F = F.*fft(U{i}(end:-1:1,:),HN,1);
                    end
                    F = ifft(F,[],1);
                    M = F(end-size(U{n},1)+1:end,:);
                else
                    % n is outside the Hankel modes
                    if n<T.dim
                        % T is not permuted and n is smaller than the Hankel
                        % modes
                        othersn = [1:n-1 n+1:T.dim-1 (T.dim+T.order):numel(T.size)];
                    else
                        % T is not permuted and n is larger than the Hankel
                        % modes
                        othersn = [1:T.dim-1 (T.dim+T.order):n-1 n+1:numel(T.size)];
                    end
                    
                    F = ones(HN,R);
                    for i = T.dim:(T.dim+T.order-1)
                        F = F.*fft(U{i},HN,1);
                    end
                    tmp = ifft(F);
                    tmp = T.val.'*tmp;
                    if ~isempty(othersn)
                        % If there are other modes besides the tensorized modes
                        if n<T.dim, contractmodes = [1:n-1 n+1:numel(T.subsize.other)];
                        else contractmodes = [1:(n-T.order-1) (n-T.order+1):numel(T.subsize.other)];
                        end
                        if n~=0, M = zeros(size(U{n}));
                        else M = zeros(1,size(U{1},2));
                        end
                        
                        for r = 1:R
                            newT = reshape(tmp(:,r),[T.subsize.other 1]);
                            ur = cellfun(@(x) x(:,r),U(othersn),'UniformOutput',false);
                            M(:,r) = contract(newT,ur,contractmodes);
                        end
                    else M = tmp;
                    end
                end
            end
            
        case 'loewner'
            if T.order>2 || ~T.isequidistant
                % Not yet implemented!
                % M = mtkrprod(ful(T),U,n,conjugate);
                % return;
                error('mtkrprod:loewner','Loewner structure is not supported yet; use the ful version instead.');
            end
            if conjugate, U = cellfun(@(x) conj(x),U,'UniformOutput',false); end
            v = T.structure.v;
            if T.ispermuted
                T.dim = 1;
            end
            
            others = [1:T.dim-1 (T.dim+T.order):numel(T.size)];
            
            f = T.val(T.ind{1},:);
            g = T.val(T.ind{2},:);
            
            if n==T.dim
                % Using Khatri-Rao product for other dimensions
                kru = kr(U{others(end:-1:1)});
                tmp = f*kru;
                
                % Using FFT for tensorized dimensions
                DB = ifft(bsxfun(@times,fft(v(end:-1:1),[],1),fft(U{T.dim+1},numel(v),1)),[],1);
                DB = DB(end-numel(v)+size(U{T.dim+1},1):end,:);
                M1 = DB.*tmp;
                
                tmp = (g*kru).*U{T.dim+1};
                M2 = ifft(bsxfun(@times,fft(v(end:-1:1),[],1),fft(tmp,numel(v),1)),[],1);
                M2 = M2(end-numel(v)+size(tmp,1):end,:);
                
                M = M1-M2;
            elseif n==T.dim+1
                % Using Khatri-Rao product for other dimensions
                kru = kr(U{others(end:-1:1)});
                
                % Using FFT for tensorized dimensions
                tmp = (f*kru).*U{T.dim};
                M1 = ifft(bsxfun(@times,fft(v,[],1),fft(tmp,numel(v),1)),[],1);
                M1 = M1(end-numel(v)+size(tmp,1):end,:);
                
                tmp = g*kru;
                DB = ifft(bsxfun(@times,fft(v,[],1),fft(U{T.dim},numel(v),1)),[],1);
                DB = DB(end-numel(v)+size(U{T.dim},1):end,:);
                M2 = DB.*tmp;
                
                M = M1-M2;
            else
                if n<T.dim
                    % T is not permuted and n is smaller than the Hankel
                    % modes
                    othersn = [1:n-1 n+1:T.dim-1 (T.dim+T.order):numel(T.size)];
                else
                    % T is not permuted and n is larger than the Hankel
                    % modes
                    othersn = [1:T.dim-1 (T.dim+T.order):n-1 n+1:numel(T.size)];
                end
                
                % Using FFT for tensorized dimensions
                DB = ifft(bsxfun(@times,fft(v(end:-1:1),[],1),fft(U{T.dim+1},numel(v),1)),[],1);
                DB = DB(end-numel(v)+size(U{T.dim+1},1):end,:);
                DBA = DB.*U{T.dim};
                
                DA = ifft(bsxfun(@times,fft(v,[],1),fft(U{T.dim},numel(v),1)),[],1);
                DA = DA(end-numel(v)+size(U{T.dim},1):end,:);
                DAB = DA.*U{T.dim+1};
                
                tmp = f.'*DBA-g.'*DAB;
                
                if ~isempty(othersn)
                    % Using Khatri-Rao product for other dimensions
                    if n<T.dim, contractmodes = [1:n-1 n+1:numel(T.subsize.other)];
                    else contractmodes = [1:(n-T.order-1) (n-T.order+1):numel(T.subsize.other)];
                    end
                    
                    if n~=0, M = zeros(size(U{n}));
                    else M = zeros(1,size(U{1},2));
                    end
                    for r = 1:R
                        newT = reshape(tmp(:,r),[T.subsize.other 1]);
                        ur = cellfun(@(x) x(:,r),U(othersn),'UniformOutput',false);
                        M(:,r) = contract(newT,ur,contractmodes);
                    end
                else M = tmp;
                end
                
            end
        otherwise
            error('mtkrprod:notImplemented', 'This is not yet implemented');
    end
    if nwas0,
        if conjugate, M = sum(M.*conj(U{1}),1);
        else M = sum(M.*U{1},1); end
    end
else
    % Determine if large-scale version of the algorithm should be executed.
    left = n == N || (n > 1 && size_tens(1) > size_tens(end));
    if left, mem = 8*R*prod(size_tens(2:end));
    else mem = 8*R*prod(size_tens(1:end-1));
    end
    largescale = mem > 2e9 || (isstruct(data) && isempty(data.matrix));
    if largescale && ~isstruct(data)
        error('mtkrprod:mem',['Intermediate result is too large. Try ' ...
            'to format data with fmt first so that the large-scale ' ...
            'version of the algorithm can be applied.']);
    end
    
    if largescale
        
        % The large scale implementation fires if the required memory would
        % otherwise be larger than 2GB and assumes the dataset has been
        % formatted by fmt.
        idx = [1:n-1 n+1:N];
        if conjugate, U = cellfun(@conj,U,'UniformOutput',false); end
        if n == 0, M = zeros(1,R);
        else M = zeros(size_tens(n),R); end
        for r = 1:R
            s = data.val;
            for j = 1:length(idx)
                s = s.*U{idx(j)}(data.sub{idx(j)},r);
            end
            if n == 0, M(r) = sum(s);
            else M(:,r) = accumarray(double(data.sub{n}),s,[size_tens(n) 1]);
            end
        end
        
    else
        
        if left
            
            % Apply tensor matrix product.
            if conjugate, M = U{1}'*reshape(T,size_tens(1),[]);
            else M = U{1}.'*reshape(T,size_tens(1),[]); end
            
            % Choose order of operations.
            perm = ones(1,N);
            l = 2; r = N;
            for i = 2:N
                if r == n || (l < n && size_tens(l) > size_tens(r))
                    perm(i) = l; l = l+1;
                else
                    perm(i) = r; r = r-1;
                end
            end
            
            % Apply structure-exploiting matricized tensor Khatri-Rao product.
            for i = 2:length(perm)
                mode = perm(i);
                if mode == n, continue; end
                if conjugate, Umode = conj(U{mode});
                else Umode = U{mode}; end
                tmp = M;
                M = alloc(R,size(M,2)/size_tens(mode));
                if mode < n
                    for j = 1:R
                        tmp2 = Umode(:,j).'* ...
                            reshape(tmp(j,:),size_tens(mode),[]);
                        M(j,:) = tmp2(:);
                    end
                elseif mode > n
                    for j = 1:R
                        tmp2 = reshape(tmp(j,:),[],size_tens(mode))*Umode(:,j);
                        M(j,:) = tmp2(:);
                    end
                end
            end
            
            % Permute output.
            M = M.';
            
        else
            
            % Apply tensor matrix product.
            if conjugate, M = reshape(T,[],size_tens(end))*conj(U{end});
            else M = reshape(T,[],size_tens(end))*U{end}; end
            
            % Choose order of operations.
            perm = N*ones(1,N);
            l = 1; r = N-1;
            for i = 2:N
                if r == n || (l < n && size_tens(l) > size_tens(r))
                    perm(i) = l; l = l+1;
                else
                    perm(i) = r; r = r-1;
                end
            end
            
            % Apply structure-exploiting matricized tensor Khatri-Rao product.
            for i = 2:length(perm)
                mode = perm(i);
                if mode == n, continue; end
                if conjugate, Umode = conj(U{mode});
                else Umode = U{mode}; end
                tmp = M;
                M = alloc(size(M,1)/size_tens(mode),R);
                if mode < n
                    for j = 1:R
                        tmp2 = Umode(:,j).'* ...
                            reshape(tmp(:,j),size_tens(mode),[]);
                        M(:,j) = tmp2(:);
                    end
                elseif mode > n
                    for j = 1:R
                        tmp2 = reshape(tmp(:,j),[],size_tens(mode))*Umode(:,j);
                        M(:,j) = tmp2(:);
                    end
                end
            end
            
        end
        
    end
    
end
end