function [H,hstruct] = hankelize(X,varargin)
%HANKELIZE Hankelization of vectors, matrices or tensors.
%   H = HANKELIZE(X) with a vector X of length N constructs a Hankel matrix
%   with X(1:ceil(N/2)) on the first column and X(ceil(N/2):end) on the
%   last row.
%
%   H = HANKELIZE(X) with a matrix X constructs a Hankel matrix for every
%   column, and stacks these matrices along the third mode of a third-order
%   tensor.
%
%   H = HANKELIZE(X) with a tensor X constructs a Hankel matrix for every
%   mode-1 fiber, and stacks them along the higher-order modes of X.
%   
%   [H,S] = HANKELIZE(X,...) returns a structure array S as an efficient
%   representation of the Hankel tensor. S can be used for efficient
%   computations such as tensor decompositions and visualizations.
%
%   H = HANKELIZE(...,'Order',K) constructs a higher-order Hankel tensor of
%   order K (instead of a Hankel matrix) for every mode-1 column/fiber of
%   X. Consider Y to be a mode-1 fiber of X. The corresponding higher-order
%   Hankel tensor will have Y(1:ceil(N/K)) in its first mode-1 fiber, and
%   Y(ceil((n-1)*N/K):ceil(n*N/K)) in the last mode-n fiber for 1<n<=K. By
%   default: K=2.
%   
%   H = HANKELIZE(...,'Dim',n) constructs a Hankel matrix or higher-order
%   Hankel tensor (depending on the order) for every mode-n fiber. By
%   default: n=1.
%
%   H = HANKELIZE(X,...,'Full',full) avoids the construction of the full
%   Hankel tensor if Full is false. H is then the efficient representation
%   of the Hankel tensor. This can be useful if the Hankel tensor is
%   large. If Full is 'auto' (by default), the full tensor is returned if
%   the storage of the tensor would not exceed the value of the FullLimit
%   option; otherwise, the efficient representation is returned. If Full is
%   true, the full tensor is always returned.
%
%   HANKELIZE(X,'key',value) or HANKELIZE(X,options) can be used to
%   pass the following options:
%   
%   - Ind:      Instead of using the segments Y(ceil((k-1)*N/K):
%               ceil(k*N/K)) of every fiber Y in X for 1<=k<=K, the
%               segments Y(1:ind(1)), Y(ind(k-1)+1:ind(k)) for 1<k<K, and
%               Y(ind(K-1)+1:end) are used. Ind is a vector of length
%               K-1.
%
%   - Sizes:    Instead of setting the breakpoints of the segments with the
%               Ind option, the user can indicate the size of the Hankel
%               tensor using the Sizes option. Sizes is a vector of length
%               K-1. The Kth element is determined as N-sum(sizes)+K-1.
%
%   - PermToFirst: The Hankel matrices or Hankel tensors are inserted at
%               the first modes, instead of at the original tensorization
%               mode. This can be useful for decompositions such as the
%               decomposition in multilinear rank-(Lr,Lr,1) terms, where Lr
%               may correspond to the rank of a Hankel representation. By
%               default: false.
%
%   - FullLimit: The storage limit to construct and return the full
%               tensor or to return the efficient representation, expressed
%               in GB. By default: 1.
%
%   See also dehankelize, loewnerize, segmentize, decimate

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] O. Debals, L. De Lathauwer, "Stochastic and Deterministic
%       Tensorization for Blind Signal Separation," Latent Variable
%       Analysis and Signal Separation, Springer Berlin / Heidelberg, Vol.
%       9237, 2015, pp. 3-13.
%
%   Version History:
%   - 2014/12/02    OD      Initial version
%   - 2015/12/14    OD      Higher-order Hankelization
%   - 2016/01/12    OD      Generalized the code to allow other modes to be
%                           Hankelized, ...
%   - 2016/02/08    OD      Changed behavior when using full=false and
%                           added storage threshold

p = inputParser();
p.addOptional('Order', 2);
p.addOptional('Dim', 1);
p.addOptional('Sizes',NaN);
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

if options.Dim>ndims(X), error('hankelize:wrongdim','The given dimension is too large!'); end

% Row to column vector
if isdefault.('Dim') && isvector(X), X = X(:); end

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

if options.Order<2, error('hankelize:smallorder','The order must be larger than or equal to 2!'); end

% Set PermToFirst
if isdefault.('PermToFirst'), options.PermToFirst = false; end
if ~islogical(options.PermToFirst), error('hankelize:permtofirst','The option PermToFirst should be a boolean!'); end
if options.PermToFirst, options.PermToFirst = 1;
else options.PermToFirst = options.Dim;
end

% Set the default values
if isdefault.('Order')
    if isnan(options.Ind)
        if isnan(options.Sizes)
            options.Ind = ceil((1:options.Order-1)*N/options.Order);
        else
            options.Order = numel(options.Sizes)+1;
            options.Ind = cumsum(options.Sizes)-(0:options.Order-2);
        end
    else
        options.Order = numel(options.Ind)+1;
        if all(~isnan(options.Sizes)) && (numel(options.Ind)~=numel(options.Sizes) || ...
                any(options.Ind~=cumsum(options.Sizes)-(0:options.Order-2)))
            error('hankelize:indsizes','The option fields Ind and Sizes do not correspond!');
        end
    end
else
    if isnan(options.Ind)
        if isnan(options.Sizes)
            options.Ind = ceil((1:options.Order-1)*N/options.Order);
        else
            if numel(options.Sizes)~=options.Order-1
                error('hankelize:sizesorder','The number of sizes should be equal to order-1');
            end
            options.Ind = cumsum(options.Sizes)-(0:options.Order-2);
        end
    else
        if all(~isnan(options.Sizes)) && (numel(options.Ind)~=numel(options.Sizes) || ...
                any(options.Ind~=cumsum(options.Sizes)-(0:options.Order-2)))
            error('hankelize:indsizes','The option fields ind and sizes do not correspond!');
        end 
    end
end
if(options.Order~=numel(options.Ind)+1)
    error('hankelize:indorder','The number of indices should be equal to order-1!');
end

% Determine the size of the Hankel matrices/tensors
size_hankel = [options.Ind N]-[0 options.Ind-1];
if any(size_hankel<0)
    error('hankelize:wrongind','The indices should increase monotonically between 0 and N!');
end
if options.Order>N
    error('hankelize:orderN','The order should not be larger than the number of elements in the tensorized dimension!');
end

% Determine the storage needed
whosX = whos('X');
totalbytes = whosX.bytes/size(X,options.Dim)*prod(size_hankel);
if ischar(options.Full) && strcmp(options.Full,'auto')
    if totalbytes > options.FullLimit*2^30
        options.Full = false;
        warning('hankelize:fullLimit',['The fullLimit has been reached, ',...
            'and instead of the dense tensor, the efficient ',...
            'representation is returned. The dense tensor can be obtained ',...
            'by setting the ''full'' option to true.']);
    else options.Full = true;
    end
end

% Determine the permutation vector applied afterwards
repermorder = [options.Order+(1:options.PermToFirst-1) 1:options.Order options.Order+options.PermToFirst:numel([size_hankel,size_other])];

if nargout>1 || ~options.Full;
    % Construct the efficient representation
    hstruct = struct;
    hstruct.type = 'hankel';
    hstruct.val = Xmat;
    hstruct.dim = options.Dim;
    hstruct.order = options.Order;
    hstruct.ind = options.Ind;
    hstruct.ispermuted = options.PermToFirst~=options.Dim;
    hstruct.repermorder = repermorder;
    hstruct.size = [size_hankel size_other ...
        ones(1,numel(hstruct.repermorder)-numel([size_hankel size_other]))];
    hstruct.size = hstruct.size(hstruct.repermorder);
    hstruct.subsize.hankel = size_hankel;
    hstruct.subsize.other = size_other;
end

if options.Full
    % Construct the (block-)Hankelized data
    allind2 = subhankel(size_hankel);
    H = Xmat(allind2,:);
    H = reshape(H,[size_hankel,size_other]);
    % Permute the data to the specific dimensions
    if options.PermToFirst~=1
        H = permute(H,repermorder);
    end
else
    % Return the efficient representation
    H = hstruct;
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