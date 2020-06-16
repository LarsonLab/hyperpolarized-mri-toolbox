function [L,lstruct] = loewnerize(X,varargin)
%LOEWNERIZE Loewnerization of vectors, matrices or tensors.
%   L = LOEWNERIZE(X) with a vector X of length N constructs a Loewner
%   matrix with
%  
%       L(i,j) = (X(p1(i))-X(p2(j))/(t(p1(i))-t(p2(j)))
%   
%   and with p1 = 1:2:N, p2 = 2:2:N and t = 1:N.
%
%   L = LOEWNERIZE(X) with a matrix X constructs a Loewner matrix for every
%   column, and stacks these matrices along the third mode of a third-order
%   tensor.
%
%   L = LOEWNERIZE(X) with a tensor X constructs a Loewner matrix for every
%   mode-1 fiber, and stacks them along the higher-order modes of X.
%   
%   [L,S] = LOEWNERIZE(X,...) returns a structure array S as an efficient
%   representation of the Loewner tensor. S can be used for efficient
%   computations such as tensor decompositions and visualizations.
%
%   L = LOEWNERIZE(...,'Order',K) constructs a higher-order Loewner tensor
%   of order K (instead of a Loewner matrix) for every mode-1 column/fiber
%   of X with elements
%       
%       L(i1,...,iK) = 
%           SUM(X(pk(ik))/PROD([(t(ik)-t(im)]],m=1,...,K,m<>k),k=1,...,K)
%   
%   with p1= 1:K:N, p2 = 2:K:N, ..., pK = K:K:N. By default: K=2.
%   
%   L = LOEWNERIZE(...,'Dim',n) constructs a Loewner matrix or higher-order
%   Loewner tensor (depending on the order) for every mode-n fiber. By
%   default: n=1.
%   
%   L = LOEWNERIZE(X,...,'full',full) avoids the construction of the full
%   Loewner tensor if Full is false. L is then the efficient representation
%   of the Loewner tensor. This can be useful if the Loewner tensor is
%   large. If Full is 'auto' (by default), the full tensor is returned if
%   the storage of the tensor would not exceed the value of the FullLimit
%   option; otherwise, the efficient representation is returned. If Full is
%   true, the full tensor is always returned.
%
%   LOEWNERIZE(X,'key',value) or LOEWNERIZE(X,options) can be used to
%   pass the following options:
%   
%   - T:        The point set t is used instead of t=1:N.
%   
%   - Ind:      Instead of using the default partitionings p1=1:K:N, ...,
%               pK=K:K:N, the partionings Ind{1}, ..., Ind{K} are used. If
%               Ind is a cell of length K-1, an additional Kth entry is
%               added. This entry contains the index values between 1 and N
%               which do not occur in Ind{1},...,Ind{K-1}.
%
%   - PermToFirst: The Löwner matrices or Löwner tensors are inserted at
%               the first modes, instead of at the original Löwnerized
%               dimension. This is useful when applying decompositions such
%               as the LL1 decomposition, for which then PermToFirst can be
%               set to true. By default: false.
%
%   - FullLimit: The limit of the storage size of the tensor determining
%               whether the (large) tensor is constructed and returned or
%               whether the efficient representation is returned. This
%               takes into account the temporary storage during the
%               tensorization step as well. Expressed in GB. By default:
%               0.5.
%
%   See also deloewnerize, hankelize, segmentize, decimate

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
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
%   - 2015/04/05    OD      Initial version
%   - 2016/01/12    OD      General higher-order Loewnerization and various
%                           options
%   - 2016/02/15    OD      Improved the memory footprint of the method.


p = inputParser();
p.addOptional('T',NaN);
p.addOptional('Order', 2);
p.addOptional('Dim', 1);
p.addOptional('Ind', NaN);
p.addOptional('PermToFirst',false);
p.addOptional('Full','auto');
p.addOptional('FullLimit',1);
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;

% Processing inputs
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);

if options.Dim>ndims(X), error('loewnerize:wrongdim','The given dimension is too large!'); end
if ~iscell(options.Ind) && ~any(isnan(options.Ind)), options.Ind = {options.Ind}; end

% Row to column vector
if isdefault.('Dim') && isvector(X), X = X(:); end

options.Ind = options.Ind(:).';
sx = size(X);
N = sx(options.Dim);
size_other = sx([1:options.Dim-1 options.Dim+1:end]);
if size_other == 1, size_other = []; end

% Matricize the data
if ismatrix(X)
    if options.Dim==2, Xmat = X.';
    else Xmat = X; end
else
    Xmat = tens2mat(X,options.Dim);
end

if iscell(options.Ind)
    if any(cell2mat(options.Ind)<1), error('loewnerize:smallind','Indices should be larger than 1'); end
    if any(cell2mat(options.Ind)>N), error('loewnerize:largeind','Indices should be smaller than the number of points in the given dimension'); end
end
if options.Order<2, error('loewnerize:smallorder','The order must be larger than or equal to 2!'); end

% Set PermToFirst
if isdefault.('PermToFirst'), options.PermToFirst = false; end
if ~islogical(options.PermToFirst), error('loewnerize:permtofirst','The option PermToFirst should be a boolean!'); end
if options.PermToFirst, options.PermToFirst = 1;
else options.PermToFirst = options.Dim;
end

% Set the default values
if isdefault.('T')
    options.T = (1:size(X,options.Dim)).';
end
t = options.T(:);
if N~=numel(t), error('loewnerize:wrongt','The number of abscissae does not agree with the number of data points'); end

if isdefault.('Ind')
    options.Ind = arrayfun(@(x) x:options.Order:N,1:options.Order,'UniformOutput',false);
else
    u = unique(cell2mat(options.Ind));
    nbind = numel(cell2mat(options.Ind));
    if numel(u)~=nbind
        error('loewnerize:overlappingind','Each index should only appear once across all of the index sets!');
    end
    if isdefault.('Order')
        if numel(u)~=N,
            options.Order = numel(options.Ind)+1;
            options.Ind{end+1} = 1:N;
            options.Ind{end}(u) = [];
        else
            options.Order = numel(options.Ind);
        end
    else
        if numel(u)~=N && options.Order~=numel(options.Ind)+1
            error('loewnerize:wrongorderind1','The order and indices do not agree!');
        elseif numel(u)==N && options.Order==numel(options.Ind)+1
            error('loewnerize:wrongorderind2','The order and indices do not agree!');
        end
        if numel(u)==N && u(1)==1 && u(N) == N
            % Correct indices
        else
            options.Ind{end+1} = 1:N;
            options.Ind{end}(u) = [];
        end
    end
end

% Determine the size of the Hankel matrices/tensors
size_loewner = cellfun(@(x) numel(x),options.Ind);

% Determine the storage needed
whosX = whos('X');
totalbytes = whosX.bytes/size(X,options.Dim)*prod(size_loewner);
if ischar(options.Full) && strcmp(options.Full,'auto')
    if totalbytes > options.FullLimit*2^30
        options.Full = false;
        warning('loewnerize:fullLimit',['The fullLimit has been reached, ',...
            'and instead of the dense tensor, the efficient ',...
            'representation is returned. The dense tensor can be obtained ',...
            'by setting the ''full'' option to true.']);
    else options.Full = true;
    end
end

% Determine the permutation vector applied afterwards
repermorder = [options.Order+(1:options.PermToFirst-1) 1:options.Order options.Order+options.PermToFirst:numel([size_loewner,size_other])];

if nargout>1 || ~options.Full
    % Construct the efficient representation
    lstruct = struct;
    lstruct.t = t;
    lstruct.dim = options.Dim;
    lstruct.ind = options.Ind;
    lstruct.order = options.Order;
    lstruct.val = Xmat;
    lstruct.type = 'loewner';
    lstruct.subsize.loewner = size_loewner;
    lstruct.subsize.other = size_other;
    lstruct.ispermuted = options.PermToFirst~=options.Dim;
    lstruct.repermorder = repermorder;
    lstruct.size = [size_loewner size_other];
    lstruct.size = [size_loewner size_other ...
        ones(1,numel(lstruct.repermorder)-numel([size_loewner size_other]))];
    lstruct.size = lstruct.size(lstruct.repermorder);
    
    % Determine if the point set T and the partitioned point sets are
    % equidistant
    if lstruct.order==2
        if all(cellfun(@(x) isequid(t(x)),lstruct.ind))
            lstruct.isequidistant = true;
            v1 = 1./(t(lstruct.ind{1}(end:-1:1))-t(lstruct.ind{2}(1)));
            v2 = 1./(t(lstruct.ind{1}(1))-t(lstruct.ind{2}(2:end)));
            lstruct.structure.v = [v1;v2];
        else
            lstruct.isequidistant = false;
        end
    else
        lstruct.isequidistant = false;
    end
end

if options.Full
    % Construct the (block-)Löwnerized data
    Xmatpart = cellfun(@(x) Xmat(x,:),options.Ind,'UniformOutput',0);
    tpart = cellfun(@(x) t(x),options.Ind,'UniformOutput',0);
    
    if isempty(size_other), size_flat_tens = [size_loewner 1];
    else size_flat_tens = [size_loewner size(Xmat,2)];
    end
    tmpind = [1 3:numel(size_flat_tens)];
    
    tpperm = cell(numel(size_loewner),1);
    for s = 1:numel(size_loewner)
        if s>1, tpperm{s} = permute(tpart{s},[2:s 1]);
        else tpperm{s} = tpart{s}; end
    end
    
    % Iteratively increase the order of the data
    for s = 1:numel(size_loewner)
        if s<numel(size_flat_tens), order = [tmpind(2:s) 1 tmpind(s+1:end) 2];
        else order = [tmpind(2:end) 1 2];
        end
        Xmpperm = permute(Xmatpart{s},order);
        size_dim = size_flat_tens; size_dim(s) = 1;
        tmp = Xmpperm;
        for r = [1:s-1 s+1:numel(size_loewner)]
            size_dimr = ones(size(size_dim)); size_dimr(r) = size_dim(r);
            tmp = repmat(tmp,size_dimr);
            tmp = bsxfun(@rdivide,tmp,bsxfun(@minus,tpperm{s},tpperm{r}));
        end
        if s==1, L = tmp; else L = L+tmp; end
    end
    L = reshape(L,[size_loewner size_other]);
    
    % Permute the data to the specific dimensions
    if options.PermToFirst~=1
        L = permute(L,repermorder);
    end
else
    % Return the efficient representation
    L = lstruct;
end

end

function f=isequid(x)
% Verify if x contains equidistant points
dx = x(2:end)-x(1:end-1);
f = all(abs(dx-dx(1))<2*eps);
end