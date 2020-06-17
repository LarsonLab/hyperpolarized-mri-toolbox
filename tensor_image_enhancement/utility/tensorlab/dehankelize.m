function X = dehankelize(H,varargin)
%DEHANKELIZE Recover the generating vector(s) from an (approximate) Hankel
%matrix/tensor.
%   X = DEHANKELIZE(H) converts an (approximate) Hankel matrix (or tensor)
%   H of size I1 x I2 x ... x IN into a column vector X of length
%   I1+I2+...+IN-N+1 by averaging along the anti-diagonal
%   (anti-hyperplanes) of the Hankel matrix (tensor).
%
%   X = DEHANKELIZE(...,'Order',K) converts an N-th order tensor of size
%   I1xI2x...xIN into a tensor of size (I1+I2+...+IK-K+1) x I_(K+1) x ... x
%   I_N by averaging each anti-diagonal or anti-hyperplane in the first K
%   modes.
%
%   X = DEHANKELIZE(...,'Dims',n) averages X along the anti-diagonals or
%   anti-hyperplanes in the modes indicated by n. The returned tensor has
%   order N-length(n)+1.
%
%   X = DEHANKELIZE(...,'Dim',d,'Order',K) is equivalent to
%   DEHANKELIZE(H,'Dims',d:d+K-1).
%
%   DEHANKELIZE(H,'key',value,...) or DEHANKELIZE(X,options) can be used to
%   pass the following options:
%
%   - Method:        Indicates the dehankelization method:
%                    - If 'fibers', the first fiber from the first mode,
%                      and the last fibers from the other modes are
%                      extracted.
%                    - If 'mean' or @mean (default), the anti-diagionals or
%                      anti-hyperplanes are averaged.
%                    - If a function handle (such as @median), this
%                      function is applied on the anti-diagonals or the
%                      anti-hyperplanes.
%                    Alternatives of 'mean' can only be used when H is a
%                    full tensor.
%
%   - PermToDim:     Permutes the result such that the PermToDim is the detensorized
%                    mode. By default: Dims(1) or Dim.
%
%   - L:             When H is a polyadic representation, L is an array
%                    determining the number of columns of the factor
%                    matrices for each dehankelization. By default, L is
%                    equal to the total number of rank-1 terms.
% 
%   See also hankelize, deloewnerize, desegmentize, dedecimate
   
%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] O. Debals, L. De Lathauwer, "Stochastic and Deterministic
%       Tensorization for Blind Signal Separation," Latent Variable
%       Analysis and Signal Separation, Springer Berlin / Heidelberg, Vol.
%       9237, 2015, pp. 3-13.
%
%   Version History:
%   - 2015/11/18   OD      Optimized version

p = inputParser();
p.addOptional('Dim',1);
p.addOptional('Dims',NaN);
p.addOptional('Order',2);
p.addOptional('Method',@mean);
p.addOptional('PermToDim',NaN);
p.addOptional('Rank',NaN);
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;

% Processing inputs
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);

if ~isdefault.('Dims')
    % dims is set
    
    if ~isdefault.('Dim')
        % dim is set
        error('dehankelize:dimanddims','Using both the ''Dim'' and ''Dims'' option arguments is invalid!');
    end
    
    if ~isdefault.('Order')
        % order is set
        if options.Dims+options.Order-1>getorder(H)
            error('dehankelize:dimsorder','The given order is not consistent with the dimensions of H!');
        end
        if numel(options.Dims)==1
            options.Dims = options.Dims:options.Dims+options.Order-1;
        end
    else
        % order is not set
        if any(options.Dims>getorder(H))
            error('dehankelize:dims','The given dehankelization dimensions are not consistent with the dimensions of H!');
        end
        options.Order = numel(options.Dims);
    end
    
    if numel(options.Dims)~=options.Order
        error('dehankelize:dimsandorder',...
            'The number of detensorized dimensions should be equal to the order!');
    end
else
    % dims is not set
    
    if ~isdefault.('Dim')
        % dim is set
        if options.Dim+options.Order-1>getorder(H)
            error('dehankelize:dim','The given dehankelization dimensions are not consistent with the dimensions of H!');
        end
    end
    
    if ~isdefault.('Order') && options.Order>getorder(H)
        error('dehankelize:order','The given order is not consistent with the dimensions of H!');
    end
    
    options.Dims = options.Dim:(options.Dim+options.Order-1);
end

if ~strcmp(getstructure(H),'cpd') && ~isdefault.('Rank')
    error('dehankelize:rank','The rank option is not supported when H is not a CPD!');
end

% Dehankelization
dims = options.Dims;
switch getstructure(H)
    case 'full'
        % Given a full tensor
        
        sh = size(H);
        size_hankel = sh(dims);
        
        if strcmp(options.Method,'fibers')
            % Recover the data by extracting single fibers
            subs = repmat({':'},1,ndims(H));
            subs(dims(2:end)) = {1};
            X = cell(options.Order,1);
            for i = 1:options.Order
                if i~=1
                    subs(dims([i-1 i])) = {size(H,dims(i-1)),2:size(H,dims(i))};
                end
                tmp = H(subs{:});
                sizes = [size(tmp) ones(1,ndims(H)-ndims(tmp))];
                permvec = 1:ndims(H); permvec(dims([1 i])) = permvec(dims([i 1]));
                tmp = permute(tmp,permvec);
                sizes = sizes(permvec);
                sizes(permvec(dims([1:i-1 i+1:numel(dims)]))) = [];
                if numel(sizes)==1, sizes(2) = 1; end
                if numel(tmp)~=0, 
                    X{i} = reshape(tmp,sizes);
                else
                    X{i} = [];
                end
            end
            signaldim = dims(1)-sum(dims(1)>dims(2:end));
            % Concatenate data
            X = cat(signaldim,X{:});
        else
            % Use another technique
            if ischar(options.Method), options.Method = str2func(options.Method); end
            Hmat = tens2mat(H,dims);
            indf = cumsum(size_hankel)-(0:1:options.Order-1);
            N = indf(end);
            
            w = hankelfreq(options.Order,size_hankel,N);
            sampleind = subhankel(size_hankel);
            
            if isequal(options.Method,@mean)
                % Extract means by performing a two-step procedure to
                % improve accuracy
                optionsquick = rmfield(options,p.UsingDefaults);
                optionsquick.PermToDim = 1;
                optionsquick.Method = 'fibers';
                xestimate = tens2mat(dehankelize(H,optionsquick),1);
                for k = 1:size(Hmat,2)
                    diff = Hmat(:,k)-xestimate(sampleind,k);
                    tmp = accumarray(sampleind,diff,[],@sum)./w;
                    xestimate(:,k) = xestimate(:,k) + tmp;
                end
            else
                % Extract something else than means (such as median)
                xestimate = nan(N,size(Hmat,2));
                for k = 1:size(Hmat,2)
                    xestimate(:,k) = accumarray(sampleind,Hmat(:,k),[],options.Method);
                end
            end
            
            sx = sh; sx(dims) = NaN;
            sx(dims(1)) = N;
            sx(isnan(sx)) = [];
            
            signaldim = dims(1)-sum(dims(1)>dims(2:end));
            X = mat2tens(xestimate,sx,signaldim);
        end
    case 'hankel'
        % Given an efficient (block-)Hankel representation
        
        if ~isdefault.('Dims') || ~isdefault.('Order')
            error('dehankelize:hankel','The dims and order are used from the Hankel field')
        end
        if ~isvector(H.val)
            X = reshape(H.val,[size(H.val,1) H.subsize.other]);
            X = ipermute(X,[H.dim 1:H.dim-1 H.dim+1:numel(H.size)]);
        else
            if H.dim==1, X = H.val;
            elseif H.dim==2, X = H.val.';
            else X = ipermute(H.val,[H.dim 1:H.dim-1]);
            end
        end
        signaldim = H.dim;
    case 'cpd'
        % Given a representation of a tensor in rank-1 terms
        
        if any(cellfun('size', H, 2) ~= size(H{1},2))
            print('For a CPD size(H{n},2) should be R for all n\n');
            return;
        end
        if ~isdefault.('Rank')
            if any(options.Rank<1)
                error('dehankelize:rank1','The different ranks should be larger than 1!')
            end
            if sum(options.Rank)~=size(H{1},2)
                error('dehankelize:ranksum','The sum of the ranks should be equal the rank of the factor matrices!');
            end
        else
            options.Rank = size(H{1},2);
        end
        rank = options.Rank;
        
        others = 1:numel(H); others(dims) = [];
        
        % Apply FFT techniques to efficiently extract the data
        size_all = cellfun('size',H,1);
        size_hankel = size_all(dims);
        N = sum(size_hankel)-options.Order+1;
        F = fft(H{dims(1)},N,1);
        for i = 2:options.Order
            F = F.*fft(H{dims(i)},N,1);
        end
        X = ifft(F,[],1);
        
        w = hankelfreq(options.Order,size_hankel,N);
        X = bsxfun(@rdivide,X,w);
        
        ranke = cumsum([1 rank]);
        i = 1:ranke(end)-1;
        j = zeros(size(i)); j(ranke(1:end-1)) = 1; j = cumsum(j);
        C = sparse(i,j,1);
         
        if ~isempty(others)
            % If other non-Hankel dimensions
            
            KR = kr(H{others(end:-1:1)}).';
            if numel(rank)~=1
                Xn = zeros(size(X,1),numel(rank),size(KR,2));
                for r = 1:numel(rank)
                    Xn(:,r,:) = X(:,ranke(r):ranke(r+1)-1)*...
                        KR(ranke(r):ranke(r+1)-1,:);
                end
                X = Xn;
            else
                X = X*KR;
            end
            
            sx = size_all; sx(dims) = NaN;
            sx(dims(1)) = N;
            if numel(rank)~=1
                sx = [sx(1:dims(1)) numel(rank) sx(dims(1)+1:end)];
            end
            sx(isnan(sx)) = [];
            
            if numel(rank)~=1
                signaldim = dims(1)-sum(dims(1)>dims(2:end));
                X = reshape(X,[size(X,1)*size(X,2) size(X,3)]);
                X = mat2tens(X,sx,[signaldim signaldim+1]);
            else
                signaldim = dims(1)-sum(dims(1)>dims(2:end));
                X = mat2tens(X,sx,signaldim);
            end
        else
            X = X*C;
            signaldim = 1;
        end
        
    case {'incomplete','sparse','tt','btd','lmlragen'}
        T = ful(H);
        X = dehankelize(T,varargin{:});
    otherwise
        error('Structure not supported!')
end

% Permute the data to specific dimensions
if ~isdefault.('PermToDim')
    options.PermToDim = double(options.PermToDim);
    signaldim = double(signaldim);
    if signaldim<options.PermToDim
        permvec = [1:signaldim-1 signaldim+1:options.PermToDim signaldim options.PermToDim+1:ndims(X)];
    else
        permvec = [1:options.PermToDim-1 signaldim options.PermToDim:signaldim-1 signaldim+1:ndims(X)];
    end
    X = permute(X,permvec);
end

end

function w = hankelfreq(order,size_hankel,N)
% Construct the frequencies of the data
if order == 2
    m = min(size_hankel(1),N-size_hankel(1)+1);
    w = [1:m-1 m*ones(1,N-2*m+2) m-1:-1:1].';
else
    w = fft(ones(size_hankel(1),1),N);
    if size(w,1)==1, w = w.'; end
    for i = 2:numel(size_hankel)
        tmp = fft(ones(size_hankel(i),1),N);
        if size(tmp,1)==1, w = w.*tmp.';
        else w = w.*tmp;
        end
    end
    w = round(ifft(w));
end
end

function allind = subhankel(sizes)
% Construct the indices of the Hankel matrix/tensor
lidx = (1:sizes(1)).';
for i = 2:numel(sizes)
    ridx = 0:sizes(i)-1;
    ridx = reshape(ridx,[ones(1,i) numel(ridx)]);
    lidx = bsxfun(@plus,lidx,ridx);
end
allind = lidx(:);
end