function X = deloewnerize(L,varargin)
%DELOEWNERIZE Recover the vector(s) from an (approximate) Löwner matrix/tensor.
%   X = DELOEWNERIZE(H) converts an (approximate) Löwner matrix (or tensor) H of
%   size I1 x I2 x ... x IN into a column vector X of length I1+I2+...+IN by
%   solving the corresponding linear Löwner system.
%
%   X = DELOEWNERIZE(...,'Order',K) converts an N-th order tensor of size
%   I1xI2x...xIN into a tensor of size (I1+I2+...+IK) x I_(K+1) x ... x
%   I_N by solving the corresponding Löwner system in the first K modes.
%
%   X = DELOEWNERIZE(...,'Dims',n) solves the corresponding Löwner system in
%   the modes indicated by n. The resulting tensor has order
%   N-length(n)+1.
%
%   X = DELOEWNERIZE(...,'Dim',d,'Order',K) is equivalent to
%   DELOEWNERIZE(H,'Dims',d:d+K-1).
%
%   DELOEWNERIZE(H,'key',value,...) or DELOEWNERIZE(X,options) can be used to
%   pass the following options:
%
%   - T:        The point set T is used instead of t=1:N.
%
%   - Ind:      Instead of using the default partitionings p1=1:K:N, ...,
%               pK=K:K:N, the partionings Ind{1}, ..., Ind{K} are used. If
%               Ind is a cell of length K-1, an additional Kth entry is
%               added. This entry contains the remaining index values
%               between 1 and N.
%
%   - PermToDim: Permutes the result such that the PermToDim is the detensorized
%               mode. By default: Dims(1) or Dim.
%
%   - L:        When the input tensor L is a polyadic representation, L is
%               an array determining the number of columns of the factor
%               matrices for each deloewnerization. By default, L is equal
%               to the total number of rank-1 terms.
%
%   See also loewnerize, dehankelize, desegmentize, dedecimate

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] O. Debals, M. Van Barel, L. De Lathauwer, "Löwner-based
%       Blind Signal Separation of Rational Functions with Applications,"
%       IEEE Transactions on Signal Processing, Vol. 64, No. 8, 2016, pp.
%       1909-1918.
%   [2] O. Debals, L. De Lathauwer, "Stochastic and Deterministic
%       Tensorization for Blind Signal Separation," Latent Variable
%       Analysis and Signal Separation, Springer Berlin / Heidelberg, Vol.
%       9237, 2015, pp. 3-13.
%   [3] K. Löwner, "Über monotone matrixfunktionen," Mathematische
%       Zeitschrift, Vol. 38, 1934, pp. 177-216.
%
%   Version History:
%   - 2014/12/02   OD      Initial version
%   - 2014/12/03   OD      Optimized version
%   - 2015/11/18   OD      Optimized version

p = inputParser();
p.addOptional('Dim',1);
p.addOptional('Dims',NaN);
p.addOptional('Order',2);
p.addOptional('PermToDim',NaN);
p.addOptional('Ind',NaN);
p.addOptional('T',NaN);
p.addOptional('Rank',NaN);
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);
options.Ind = options.Ind(:).';

% Error checks
if ~isdefault.('Dims')
    % dims is set
    if ~isdefault.('Dim')
        % dim is set
        error('deloewnerize:dimanddims','Using both the ''Dim'' and ''Dims'' option arguments is invalid!');
    end
    
    if ~isdefault.('Order')
        % order is set
        if options.Dims+options.Order-1>getorder(L)
            error('deloewnerize:dimsorder','The given order is not consistent with the dimensions of L!');
        end
        if numel(options.Dims)==1
            options.Dims = options.Dims:options.Dims+options.Order-1;
        end
    else
        % order is not set
        if any(options.Dims>getorder(L))
            error('deloewnerize:dims','The given deloewnerization dimensions are not consistent with the dimensions of L!');
        end
        options.Order = numel(options.Dims);
    end
    
    if numel(options.Dims)~=options.Order
        error('deloewnerize:dimsandorder',...
            'The number of detensorized dimensions should be equal to the order!');
    end
else
    % dims is not set
    if ~isdefault.('Dim')
        % dim is set
        if options.Dim+options.Order-1>getorder(L)
            error('deloewnerize:dim','The given deloewnerization dimensions are not consistent with the dimensions of H!');
        end
    end
    
    if ~isdefault.('Order') && options.Order>getorder(L)
        error('deloewnerize:order','The given order is not consistent with the dimensions of L!');
    end
    
    options.Dims = options.Dim:(options.Dim+options.Order-1);
    
end

switch getstructure(L)
    case 'full'
        % Given a full tensor
        
        sh = size(L);
        size_loewner = sh(options.Dims);
        N = sum(size_loewner);
        
        if isdefault.('T')
            options.T = 1:N;
        end
        
        if isdefault.('Ind')
            options.Ind = arrayfun(@(x) x:options.Order:N,1:options.Order,'UniformOutput',false);
        else
            if ~iscell(options.Ind), options.Ind = {options.Ind}; end
            u = unique(cell2mat(options.Ind));
            if numel(u)~=N && options.Order~=numel(options.Ind)+1
                error('deloewnerize:wrongorderind1','The order and indices do not agree!');
            elseif numel(u)==N && options.Order==numel(options.Ind)+1
                error('deloewnerize:wrongorderind2','The order and indices do not agree!');
            end
            if numel(u)==N && u(1)==1 && u(N) == N
                % Correct indices
            else
                options.Ind{end+1} = 1:N;
                options.Ind{end}(u) = [];
            end
        end
        sizeind = cellfun('length',options.Ind);
        for i=1:numel(options.Ind)
            if sizeind(i)~=size_loewner(i)
                error('deloewnerize:indloewnersize',...
                    'The number of indices in mode %d does not match the corresponding dimension of the Löwner tensor!',i)
            end
        end
        Lmat = tens2mat(L,options.Dims);
        
        % Construct the left-hand matrix of the system
        F = struct;
        F.size = [size_loewner N];
        F.sparse = 1;
        F.sub = cell(1,options.Order+1);
        F.val = [];
        t = options.T(:);
        C = num2cell(ones(1,options.Order));
        for o = 1:options.Order
            v = [1:o-1 o+1:options.Order];
            for mi = 1:options.Order-1
                m = v(mi);
                tmp = bsxfun(@minus,t(options.Ind{o}),t(options.Ind{m}).');
                tmp = reshape(tmp,[size(tmp,1) ones(1,mi-1) size(tmp,2)]);
                C{o} = bsxfun(@times,C{o},tmp);
            end
            C{o} = 1./C{o};
            ndgridind = arrayfun(@(x) 1:x,size_loewner,'UniformOutput',false);
            
            Fsub = cell(1,options.Order+1);
            [Fsub{1:end-1}] = ndgrid(ndgridind{:});
            Fsub{end} = options.Ind{o}(Fsub{o}(:));
            F.sub = cellfun(@(x,y) [x;y(:)],F.sub,Fsub,'UniformOutput',false);
            Cm = permute(C{o},[2:o 1 o+1:getorder(C{o})]);
            F.val = [F.val;Cm(:)];
        end
        
        sx = sh; sx(options.Dims) = NaN;
        sx(options.Dims(1)) = N;
        sx(isnan(sx)) = [];
        signaldim = options.Dims(1)-sum(options.Dims(1)>options.Dims(2:end));
        G = tens2mat(F,[],getorder(F));
        
        % Solve the system, and pad with zeros
        X = G(:,1:end-options.Order+1)\Lmat;
        X = [X;zeros(options.Order-1,size(X,2))];
        
        % Normalize to obtain zero mean, ...
        if options.Order==2
            X = bsxfun(@minus,X,mean(X,1));
        else
            V = bsxfun(@power,linspace(-1,1,size(X,1)).',options.Order-2:-1:0);
            [Q,R] = qr(V,0);
            c = zeros(options.Order-1,size(X,2));
            for n = 1:size(X,2)
                c(:,n) = (R\(Q.'*X(:,n))).';
            end
            X = X-V*c;
        end
        X = mat2tens(X,sx,signaldim);
        
    case 'loewner'
        % Given an efficient (block-)Löwner representation
        
        if ~isdefault.('Dims') || ~isdefault.('Order') || ~isdefault.('Ind') || ~isdefault.('T')
            error('deloewnerize:loewner','The dims, order, ind and t are used from the Löwner field')
        end
        if ~isvector(L.val)
            X = reshape(L.val,[size(L.val,1) L.subsize.other]);
            newrepermorder = L.repermorder;
            newrepermorder(newrepermorder > 1 & newrepermorder <= L.order) = [];
            newrepermorder(newrepermorder > L.order) = newrepermorder(newrepermorder>L.order)-L.order+1;
            X = permute(X,newrepermorder);
        else
            if L.dim==1, X = L.val;
            elseif L.dim==2, X = L.val.';
            else X = permute(L.val,[L.dim 1:L.dim]);
            end
        end
        if L.ispermuted, signaldim = 1;
        else signaldim = L.dim;
        end
    case 'cpd'
        % Given a representation of a tensor in rank-1 terms
        
        if any(cellfun('size', L, 2) ~= size(L{1},2))
            print('For a CPD size(H{n},2) should be R for all n\n');
            return;
        end
        if ~isdefault.('Rank')
            if any(options.Rank<1)
                error('deloewnerize:rank1','The different ranks should be larger than 1!')
            end
            if sum(options.Rank)~=size(L{1},2)
                error('deloewnerize:ranksum','The sum of the ranks should be equal the rank of the factor matrices!');
            end
        else
            options.Rank = size(L{1},2);
        end
        rank = options.Rank;
        ranke = cumsum([1 rank]);
        if numel(rank)~=1
            Xt = cell(1,numel(rank));
            for r = 1:numel(rank)
                T = ful(cellfun(@(x) x(:,ranke(r):ranke(r+1)-1) ,L,'UniformOutput',false));
                Xt{r} = deloewnerize(T,varargin{:});
            end
            X = cat(ndims(Xt{1}),Xt{:});
            X = permute(X,[1:options.Dims(1) ndims(X) options.Dims(1)+1:ndims(X)-1]);
        else
            T = ful(L);
            X = deloewnerize(T,varargin{:});
        end
    case {'incomplete','sparse','tt','btd','lmlragen'}
        T = ful(L);
        X = deloewnerize(T,varargin{:});
    otherwise
        error('Structure not supported!')
end

% Permute the data to specific dimensions
if ~isdefault.('PermToDim')
    if signaldim<options.PermToDim
        permvec = [1:signaldim-1 signaldim+1:options.PermToDim signaldim options.PermToDim+1:ndims(X)];
    else
        permvec = [1:options.PermToDim-1 signaldim options.PermToDim:signaldim-1 signaldim+1:ndims(X)];
    end
    X = permute(X,permvec);
end
