function f = inprod(T1,T2)
%INPROD Inner product of two tensors
%   INPROD(T1,T2) computes the inner product T2(:)'*T1(:) between two
%   tensors. The tensors can structured. Currently, combinations of the
%   following types are supported: full, incomplete, sparse, cpd, lmlra, btd,
%   tt, hankel and loewner. In the case of incomplete tensors, the weighted
%   inner product is computed, in which the binary weight tensor is one for
%   known elements. In the case of two incomplete tensors, the binary weight
%   tensor is one for entries known in both tensors.
%
%   See also frob, getstructure.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2015/12/23   OD      Hankel, Loewner
% - 2015/12/11   NV      CPD, LMLRA, BTD, TT

f = nan;
type1 = getstructure(T1);
type2 = getstructure(T2);

types = {'full','incomplete','sparse','cpd', 'lmlra', 'btd', 'tt','hankel','loewner'};
i = find(cellfun(@(s) strcmp(type1, s), types));
j = find(cellfun(@(s) strcmp(type2, s), types));

if isempty(i)
    error('inprod:unknownTypes', 'The type ''%s'' is not supported', type1);
end
if isempty(j)
    error('inprod:unknownTypes', 'The type ''%s'' is not supported', type2);
end

% Check dimensions
sz1orig = getsize(T1);
sz2orig = getsize(T2);
sz1 = sz1orig(1:find(sz1orig > 1, 1, 'last'));
sz2 = sz2orig(1:find(sz2orig > 1, 1, 'last'));
if length(sz1) ~= length(sz2)
    error('inprod:dimensionMismatch', ['The orders of the tensors should ' ...
        'match: getorder(T1) ~= getorder(T2)']);
end
if any(sz1 ~= sz2)
    error('inprod:dimensionMismatch', ['The dimensions of the tensors should ' ...
        'match: getsize(T1) ~= getorder(T2)']);
end
if find(sz1orig > 1, 1, 'last') < length(sz1orig),
    T1 = simplify(T1,sz1orig,type1);
end
if find(sz2orig > 1, 1, 'last') < length(sz2orig),
    T2 = simplify(T2,sz2orig,type2);
end

% The inner product is conjugated symmetric, so only half of the combinations
% is needed.
if j > i
    f = conj(inprod(T2,T1));
    return;
end

switch type2
    case 'full'
        N2 = ndims(T2);
        switch type1
            case 'full'
                f = T2(:)'*T1(:);
            case 'incomplete'
                T1 = fmt(T1);
                f = T2(T1.ind)'*T1.val(:);
            case 'sparse'
                T2(T1.ind) = T2(T1.ind).*conj(T1.val);
                f = sum(conj(T2(:)));
            case 'cpd'
                tmp = ones(1, size(T1{1},2));
                for r = 1:size(T1{1},2)
                    u = cellfun(@(u) conj(u(:,r)), T1(1:N2), 'UniformOutput', false);
                    tmp(r) = tmp(r)*conj(contract(T2, u, 1:N2));
                end
                f = sum(tmp);
            case 'lmlra'
                S = tmprod(T2,T1{1},1:length(T1{1}),'H');
                f = conj(T1{2}(:)'*S(:));
            case 'btd'
                f = 0;
                for r = 1:length(T1);
                    f = f + inprod([{T1{r}(1:end-1)} T1{r}{end}], T2);
                end
            case 'tt'
                ttr = [1 cellfun('size', T1(2:end), 1) 1];
                N1 = min(length(T1),N2);
                size_res = sz2;
                for n = 1:N1
                    size_res(1) = ttr(n)*sz2(n);
                    T2 = reshape(T2,size_res(1),prod(size_res(2:end)));
                    T2 = reshape(T1{n},ttr(n)*sz2(n),ttr(n+1))'*T2;
                    size_res(1) = [];
                end
                f = conj(T2);
            case 'hankel'
                
                T1f = ful(T1);
                f = inprod(T1f,T2);
                
                return;
                
                % 'efficient' implementation is slower than using ful.
                
                if ~T1.ispermuted
                    T2h = dehankelize(T2,'order',T1.order,'dim',T1.dim);
                    T2h = tens2mat(T2h,[T1.dim 1:T1.dim-1 T1.dim+1:ndims(T2h)],[]);
                else
                    T2h = dehankelize(T2,'order',T1.order,'dim',1);
                    T2h = T2h(:);
                end
                size_hankel = T1.subsize.hankel;
                HN = size(T1.val,1);
                w = fft(ones(size_hankel(1),1),HN);
                if size(w,1)==1, w = w.'; end
                for i = 2:numel(size_hankel)
                    tmp = fft(ones(size_hankel(i),1),HN);
                    if size(tmp,1)==1, w = w.*tmp.';
                    else w = w.*tmp;
                    end
                end
                w = round(ifft(w));
                tmp = bsxfun(@times,T1.val,w);
                f = T2h'*tmp(:);
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
    case 'incomplete'
        T2 = fmt(T1);
        switch type1
            case {'incomplete', 'sparse'}
                T1 = fmt(T1);
                % faster setdiff
                [i,j] = sort([T1.ind(:); T2.ind(:)]);
                k = diff(i) == 0;
                rem1 = j(k);
                rem2 = j([false; k]) - length(T1.ind);
                f = T2.val(rem2)'*T1.val(rem1);
            case 'cpd'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'lmlra'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'btd'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'tt'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'hankel'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
    case 'sparse'
        T2 = fmt(T2);
        switch type1
            case 'sparse'
                T1 = fmt(T1);
                % faster setdiff
                [i,j] = sort([T1.ind(:); T2.ind(:)]);
                k = diff(i) == 0;
                rem1 = j(k);
                rem2 = j([false; k]) - length(T1.ind);
                f = T2.val(rem2)'*T1.val(rem1);
            case 'cpd'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'lmlra'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'btd'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'tt'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'hankel'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
    case 'cpd'
        R  = size(T2{1},2);  % rank of cpd
        N2 = length(T2);     % number of dimensions of T2
        
        switch type1
            case 'cpd'
                W = cellfun(@(u,v) v'*u, T1, T2, 'UniformOutput', false);
                W = prod(cat(3, W{:}),3);
                f = sum(W(:));
            case 'lmlra'
                W = cellfun(@(u,v) v'*u, T1{1}, T2, 'UniformOutput', false);
                S = T1{2};
                if R^ndims(S) > 1e8 % large scale version
                    W = cellfun(@transpose, W, 'UniformOutput', false);
                    tmp = nan(1, R);
                    for r = 1:R
                        u = cellfun(@(u) u(:,r), W, 'UniformOutput', false);
                        tmp(r) = contract(S,u(1:ndims(S)),1:ndims(S));
                        % in the case ndims(S)<ndims(T1{2})
                        tmp(r) = tmp(r)*prod(cat(1,u{ndims(S)+1:end}));
                    end
                else % small scale version
                    tmp = tmprod(S, W(1:ndims(S)), 1:ndims(S));
                    tmp = tmp(linspace(1,numel(tmp),R));
                    if ndims(S) < length(W)
                        tmp = tmp.*prod(cat(2, W{ndims(S)+1:end}),2).';
                    end
                end
                f = sum(tmp(:));
            case 'btd'
                f = 0;
                for r = 1:length(T1)
                    f = f + inprod([{T1{r}(1:end-1)}, T1{r}{end}], T2);
                end
            case 'tt'
                ttr = [1 cellfun('size', T1(2:end), 1) 1];
                size_tens = cellfun('size', T2, 1);
                N1 = min(length(T1), N2);
                tmp = cell(1, N1);
                tmp{1} = reshape(T1{1}.'*conj(T2{1}), 1, ttr(2), R);
                for n = 2:N1
                    tmp{n} = reshape(permute(T1{n}, [1 3 2]), ttr(n)*ttr(n+1), size_tens(n));
                    tmp{n} = reshape(tmp{n}*conj(T2{n}), ttr(n), ttr(n+1), R);
                end
                t = zeros(1, R);
                for r = 1:R
                    res = tmp{1}(:,:,r);
                    for n = 2:length(tmp)
                        res = res * tmp{n}(:,:,r);
                    end
                    t(r) = res;
                end
                f = sum(t(:));
            case 'hankel'
                HN = size(T1.val,1);
                
                % The tensorized dimensions
                [~,dims] = ismember(1:T1.order,T1.repermorder);
                % The non-tensorized dimensions
                others = 1:getorder(T1); others(dims) = [];
                
                T2 = cellfun(@conj,T2,'UniformOutput',false);
                % Taking HN-point FFT of the factor matrices in the
                % tensorized dimensions, and multiplying it elementwise
                % along the first dimension
                F = fft(T2{dims(1)},HN,1);
                for i = 2:numel(dims)
                    F = F.*fft(T2{dims(i)},HN,1);
                end
                % Taking IFFT
                F = ifft(F,[],1);
                
                if ~isempty(others)
                    % If there are other dimensions, we can apply
                    tmp = cpdgen([F,T2(others)]);
                else
                    % Else we sum the columns (faster than sum(...,2))
                    tmp = F*ones(size(F,2),1);
                end
                f = tmp(:).'*T1.val(:);
                
            case 'loewner'
                if T1.order>2 || ~T1.isequidistant
                    % The efficient representation is not supported for
                    % higher-order Loewner tensors, or when the abscissae
                    % are not equidistant.
                    T1f = ful(T1);
                    f = inprod(T1f,T2);
                    return
                end
                
                T2 = cellfun(@conj,T2,'UniformOutput',false);
                
                % The tensorized dimensions
                [~,dims] = ismember(1:T1.order,T1.repermorder);
                % The non-tensorized dimensions
                others = 1:getorder(T1); others(dims) = [];
                
                % Partition data in f (ind{1}) and g (ind{2})
                f = T1.val(T1.ind{1},:);
                g = T1.val(T1.ind{2},:);
                v = T1.structure.v;
                
                % Calculate the product. Basically, the product
                % [(fi-gj)/(xi-yj)]*tij is separated in fi*tij/(xi-yj) and
                % -gj*tij/(xi-yj). The division by (xi-yj) can be
                % efficiently implemented as the Cauchy matrix 1/(xi-yj) is
                % a Toeplitz matrix for equidistant points. Hence, a
                % convolution appears, and the FFT/IFFT can be used. DBA
                % can be seen as t:j/(xi-yj) and DAB can be seen as
                % ti:/(xi-y:).
                
                DB = ifft(bsxfun(@times,...
                    fft(v(end:-1:1),[],1),...
                    fft(T2{dims(2)},numel(v),1)),...
                    [],1);
                DB = DB((end-numel(v)+size(T2{dims(2)},1)):end,:);
                DBA = DB.*T2{dims(1)};
                
                DA = ifft(bsxfun(@times,...
                    fft(v,[],1),...
                    fft(T2{dims(1)},numel(v),1)),[],1);
                DA = DA(end-numel(v)+size(T2{dims(1)},1):end,:);
                DAB = DA.*T2{dims(2)};
                
                % Calculate the substraction
                sub = f.'*DBA-g.'*DAB;
                
                if ~isempty(others)
                    % If there are other dimensions, we can apply
                    tmp = kr(T2{others(end:-1:1)});
                    f = sub(:).'*tmp(:);
                else
                    % Else we sum the columns (faster than sum(...,2))
                    f = sum(sub(:));
                end
        end
    case 'lmlra'
        N2 = length(T2{1});
        switch type1
            case 'lmlra'
                W = cellfun(@(u,v) v'*u, T1{1}, T2{1}, 'UniformOutput', false);
                S = tmprod(T2{2}, W, 1:length(W), 'H');
                f = S(:)'*T1{2}(:);
            case 'btd'
                f = 0;
                for r = 1:length(T1)
                    f = f + inprod([{T1{r}(1:end-1)}, T1{r}{end}], T2);
                end
            case 'tt'
                ttr = [1 cellfun('size', T1(2:end), 1) 1];
                size_tens = cellfun('size', T2{1}, 1);
                size_core = [size(T2{2}), ones(1,N2-ndims(T2{2}))];
                N1 = min(length(T1), N2);
                tmp = cell(1, N1);
                tmp{1} = T2{1}{1}'*T1{1};
                for n = 2:N1
                    tmp{n} = reshape(permute(T1{n}, [1 3 2]), ttr(n)*ttr(n+1), size_tens(n));
                    tmp{n} = reshape(tmp{n}*conj(T2{1}{n}), ttr(n), ttr(n+1), size_core(n));
                    tmp{n} = permute(tmp{n}, [1 3 2]);
                end
                t = ttgen(tmp);
                f = T2{2}(:)'*t(:);
            case 'hankel'
                
                HN = size(T1.val,1);
                % The tensorized dimensions
                [~,idx] = sort(T1.repermorder);
                dims = idx(1:T1.order);
                % The non-tensorized dimensions
                others = 1:getorder(T1); others(dims) = [];
                
                T2{1} = cellfun(@conj,T2{1},'UniformOutput',false);
                T2{2} = conj(T2{2});
                
                % Taking HN-point FFT of the factor matrices in the
                % tensorized dimensions, and multiplying it elementwise
                % along the first dimension
                F = fft(T2{1}{dims(end)},HN,1).';
                for i = T1.order-1:-1:1
                    F = kr(F,fft(T2{1}{dims(i)},HN,1).');
                end
                % Taking IFFT
                F = F.';
                F = ifft(F,[],1);
                
                % Reshape core of T2 to the right sizes
                Ss = desegmentize(T2{2},'Dims',dims(dims<=getorder(T2{2})),'PermToDim',1);
                
                if ~isempty(others)
                    % If there are other dimensions
                    F = tmprod(Ss,{F,T2{1}{others}},1:numel(others)+1);
                    tmp = F(:);
                else
                    % If there are no other dimensions
                    tmp = F*Ss;
                end
                f = tmp.'*T1.val(:);
                
            case 'loewner'
                if T1.order>2 || ~T1.isequidistant
                    % The efficient representation is not supported for
                    % higher-order Loewner tensors, or when the abscissae
                    % are not equidistant.
                    T1f = ful(T1);
                    f = inprod(T1f,T2);
                    return
                end
                
                T2{1} = cellfun(@conj,T2{1},'UniformOutput',false);
                T2{2} = conj(T2{2});
                
                % The tensorized dimensions
                [~,dims] = ismember(1:T1.order,T1.repermorder);
                % The non-tensorized dimensions
                others = 1:getorder(T1); others(dims) = [];
                
                % Partition data in f (ind{1}) and g (ind{2})
                f = T1.val(T1.ind{1},:);
                g = T1.val(T1.ind{2},:);
                v = T1.structure.v;
                
                % Calculate the product. Basically, the product
                % [(fi-gj)/(xi-yj)]*tij is separated in fi*tij/(xi-yj) and
                % -gj*tij/(xi-yj). The division by (xi-yj) can be
                % efficiently implemented as the Cauchy matrix 1/(xi-yj) is
                % a Toeplitz matrix for equidistant points. Hence, a
                % convolution appears, and the FFT/IFFT can be used. DBA
                % can be seen as t:j/(xi-yj) and DAB can be seen as
                % ti:/(xi-y:).
                
                DB = ifft(bsxfun(@times,...
                    fft(v(end:-1:1),[],1),...
                    fft(T2{1}{dims(2)},numel(v),1)),...
                    [],1);
                DB = DB((end-numel(v)+size(T2{1}{dims(2)},1)):end,:); % I x v
                DBA = kr(DB.',T2{1}{dims(1)}.');
                DA = ifft(bsxfun(@times,...
                    fft(v,[],1),...
                    fft(T2{1}{dims(1)},numel(v),1)),[],1);
                DA = DA(end-numel(v)+size(T2{1}{dims(1)},1):end,:); % J x u
                DAB = kr(T2{1}{dims(2)}.',DA.');
                
                % Calculate the substraction
                sub = DBA*f-DAB*g;
                
                % Reshape core of T2 to the right sizes
                % Ss = desegmentize(T2{2},'dims',dims(dims<=getorder(T2{2})),'permdim',1);
                
                Ss = tens2mat(T2{2},dims(dims<=getorder(T2{2})),[]);
                
                if ~isempty(others)
                    % If there are other dimensions
%                     tmp = mtkronprod(reshape(sub,[size(sub,1) arrayfun(@(x) size(T2{1}{x},1),others)]),T2{1}(others),0);
                    tmp = sub*kron(T2{1}{others(end:-1:1)});
                    f = Ss(:).'*tmp(:);
                else
                    % If there are no other dimensions
                    f = sub(:).'*Ss(:);
                end
        end
    case 'btd'
        switch type1
            case 'btd'
                f = 0;
                for r = 1:length(T1)
                    for s = 1:length(T2)
                        f = f + inprod([{T1{r}(1:end-1)}, T1{r}{end}], ...
                            [{T2{s}(1:end-1)}, T2{s}{end}]);
                    end
                end
            case 'tt'
                f = 0;
                for r = 1:length(T2)
                    f = f + inprod(T1, [{T2{r}(1:end-1)}, T2{r}{end}]);
                end
            case 'hankel'
                f = 0;
                for r = 1:length(T2)
                    f = f + inprod(T1, [{T2{r}(1:end-1)}, T2{r}{end}]);
                end
            case 'loewner'
                f = 0;
                for r = 1:length(T2)
                    f = f + inprod(T1, [{T2{r}(1:end-1)}, T2{r}{end}]);
                end
        end
    case 'tt'
        N2 = length(T2);
        switch type1
            case 'tt'
                size_tens = cellfun('size', T1, 2);
                ttr1 = [1 cellfun('size', T1(2:end), 1) 1];
                ttr2 = [1 cellfun('size', T2(2:end), 1) 1];
                N1 = min(length(T1),N2);
                tmp = T2{1}'*T1{1};
                f = tmp(:).';
                for n = 2:N1-1
                    tmp1 = reshape(permute(T1{n},[2 1 3]),size_tens(n),ttr1(n)*ttr1(n+1));
                    tmp2 = reshape(permute(T2{n},[2 1 3]),size_tens(n),ttr2(n)*ttr2(n+1));
                    tmp = tmp2'*tmp1;
                    tmp = reshape(tmp, ttr2(n), ttr2(n+1), ttr1(n), ttr1(n+1));
                    tmp = permute(tmp, [1 3 2 4]);
                    tmp = reshape(tmp, ttr2(n)*ttr1(n), ttr2(n+1)*ttr1(n+1));
                    f = f * tmp;
                end
                tmp = conj(T2{N2}*T1{N1}');
                f = f*tmp(:);
            case 'hankel'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
    case 'hankel'
        switch type1
            case 'hankel'
                T1f = ful(T1);
                f = inprod(T1f,T2);
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
    case 'loewner'
        switch type1
            case 'loewner'
                T1f = ful(T1);
                f = inprod(T1f,T2);
        end
        
    otherwise
        error('Not implemented');
end

end

function T = simplify(T, sz, type)
n = find(sz > 1, 1, 'last');
switch type
    case 'cpd'
        T{n} = bsxfun(@times, T{n}, prod(cat(1, T{n+1:end}),1));
        T(n+1:end) = [];
    case 'lmlra'
        T{2} = tmprod(T{2}, T{1}(n+1:end), n+1:length(T{1}));
        T{1}(n+1:end) = [];
    case 'btd'
        for r = 1:length(T)
            n = find(cellfun('size',T{r}(1:end-1),1) > 1, 1, 'last');
            if n ~= length(T{r})-1
                T{r}{end} = tmprod(T{r}{end}, T{r}(n+1:end-1), n+1:length(T{r})-1);
                T{r}(n+1:end-1) = [];
            end
        end
    case 'tt'
        for k = length(T):-1:n+1
            tmp = T{k-1};
            tmp = reshape(tmp, size(T{k-1},1)*size(T{k-1},2),size(T{k-1},3));
            T{k-1} = reshape(tmp*T{k}, size(T{k-1},1),size(T{k-1},2));
        end
        T(n+1:end) = [];
end
end
