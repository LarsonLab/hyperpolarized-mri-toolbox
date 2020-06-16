function [S,sstruct] = segmentize(X,varargin)
%SEGMENTIZE Segmentation of vectors, matrices or tensors.
%   S = SEGMENTIZE(X) with a vector X of length N=IJ constructs a matrix of
%   size IxJ with the J consecutive segments of length I in its columns.
%
%   S = SEGMENTIZE(X) with a matrix X constructs a matrix from the different
%   segments of every column, and stacks these matrices along the third
%   mode of a third-order tensor.
%
%   S = SEGMENTIZE(X) with a tensor X constructs a matrix from the different
%   segments of every mode-1 fiber, and stacks them along the higher-order
%   modes of X.
%   
%   [S,str] = SEGMENTIZE(X,...) returns a structure array str as an
%   efficient representation of the output tensor. str can be used for
%   efficient computations such as tensor decompositions and
%   visualizations.
%
%   S = SEGMENTIZE(...,'Order',K) constructs a tensor of order K from each
%   fiber by reshaping the data into a tensor of size I_1x...xI_K. The
%   values I_1,...,I_K are chosen such that they are as close together as
%   possible and such that as much data of X is used as possible. By default:
%   K=2.
%
%   S = SEGMENTIZE(...,'Dim',n) constructs a matrix or higher-order tensor
%   (depending on the order) for every mode-n fiber. By default: n=1.
%
%   S = SEGMENTIZE(X,...,'Full',full) avoids the construction of the full
%   tensor if Full is false. S is then the efficient representation
%   of the output tensor. This can be useful if the output tensor is
%   large. If Full is 'auto' (by default), the full tensor is returned if
%   the storage of the tensor would not exceed the value of the FullLimit
%   option; otherwise, the efficient representation is returned. If Full is
%   true, the full tensor is always returned.
%
%   SEGMENTIZE(X,'key',value) or SEGMENTIZE(X,options) can be used to
%   pass the following options:
%
%   - Segsize:          Instead of the default values I_1,...,I_K-1, the
%                       values from Segsize are used. I_K is chosen such
%                       that as much data of X is used as possible. Segsize
%                       is a vector of length K-1. If the Nsegments option
%                       is used (see further), an error is thrown.
%
%   - Nsegments:        Instead of the default values I_2,...,I_K, the
%                       values from Nsegments are used. I_1 is chosen such
%                       that as much data of X is used as possible.
%                       Nsegments is a vector of length K-1. If the Segsize
%                       option is used as well, an error is thrown.
%
%   - Shift:            Shift is a vector of length K-1, indicating the
%                       shift between consecutive segments. By default, the
%                       elements of Shift are equal to the frame lengths
%                       such that there is no overlapping. If a shift equal
%                       to 1 is used, a Hankel structure is obtained in
%                       that mode. If all Shifts are equal to 1, a
%                       higher-order Hankel tensor is obtained from each
%                       tensorized fiber.
%
%   - PermToFirst:      The matrices or tensors are inserted at
%                       the first modes, instead of at the original
%                       tensorized mode. This can be useful for
%                       decompositions such as the decomposition in
%                       multilinear rank-(Lr,Lr,1) terms, where Lr may
%                       correspond to the rank of a Hankel representation.
%                       By default: false.
%
%   - UseAllSamples:    If true, an error is thrown if not all of the
%                       samples of X are used. By default: false.
%
%   - FullLimit:        The storage limit to construct and return the full
%                       tensor or to return the efficient representation,
%                       expressed in GB. By default: 1.
%
%   See also desegmentize, hankelize, loewnerize, decimate

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] M. Boussé, O. Debals, L. De Lathauwer, "A Tensor-Based Method for
%       Large-Scale Blind Source Separation using Segmentation," Internal
%       Report 15-59, ESAT-STADIUS, KU Leuven, Belgium, 2015.
%   [2] M. Boussé, O. Debals, L. De Lathauwer, "A novel deterministic
%       method for large-scale blind source separation," Proceedings of the
%       23rd European Signal Processing Conference (EUSIPCO, Nice, France),
%       2015
%   [3] O. Debals, L. De Lathauwer, "Stochastic and Deterministic
%       Tensorization for Blind Signal Separation," Latent Variable
%       Analysis and Signal Separation, Springer Berlin / Heidelberg, Vol.
%       9237, 2015, pp. 3-13.
%
%   Version History:
%   - 2015/01/10    OD      Initial version

p = inputParser();
p.addOptional('Order', 2);
p.addOptional('Dim', 1);
p.addOptional('Segsize', NaN);
p.addOptional('Nsegments',NaN);
p.addOptional('PermToFirst',false);
p.addOptional('Full','auto');
p.addOptional('Shift',NaN);
p.addOptional('UseAllSamples',false);
p.addOptional('FullLimit',1);
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;

% Processing inputs
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);

if options.Dim>ndims(X), error('segmentize:wrongdim','The given dimension is too large!'); end
if ~any(isnan(options.Shift)) && (any(options.Shift<0) || any(options.Shift==0)), error('segmentize:shift','The shifts cannot be negative or zero!'); end

% Row to column vector
if any(strcmp('Dim',p.UsingDefaults)) && isvector(X), X = X(:); end

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

if ~isdefault.('Segsize') && ~isdefault.('Nsegments')
    error('segmentize:segsizeandNsegments','Only one of the options Segsize and Nsegments is allowed to be used at the same time.');
end

% Set the default values
if any(strcmp('Order',p.UsingDefaults))
    % Default order
    if isnan(options.Segsize)
        % Default segment size
        if isnan(options.Nsegments)
            % Default number of segments
            if isnan(options.Shift)
                options.Order = 2;
            else
                options.Order = numel(options.Shift)+1;
            end
            options.Segsize = ccd(N,[],options.Order);
        else
            % Number of segments set
            options.Order = numel(options.Nsegments)+1;
        end
    else
        % Segment size set
        options.Order = numel(options.Segsize)+1;
    end
else
    % Order set
    if isnan(options.Segsize)
        % Default segment size
        if isnan(options.Nsegments)
            % Default number of segments
            options.Segsize = ccd(N,[],options.Order);
        elseif(options.Order~=numel(options.Nsegments)+1)
            error('segmentize:wrongorderNsegments','The order and number of segment numbers do not correspond!');
        end
    elseif(options.Order~=numel(options.Segsize)+1)
        error('segmentize:wrongordersegsize','The order and number of segment sizes do not correspond!');
    end
end

if options.Order<2, error('segmentize:smallorder','The order must be larger then or equal to 2!'); end

% Set PermToFirst
if isdefault.('PermToFirst'), options.PermToFirst = false; end
if ~islogical(options.PermToFirst), error('segmentize:permtofirst','The option PermToFirst should be a boolean!'); end
if options.PermToFirst, options.PermToFirst = 1;
else options.PermToFirst = options.Dim;
end

if all(~isnan(options.Shift)) && (numel(options.Shift)~=options.Order-1) && ~any(strcmp('Order',p.UsingDefaults))
    error('segmentize:shiftorder','The number of shifts is not equal to the order minus one!');
end
if any(isnan(options.Segsize))
    if all(~isnan(options.Nsegments)) && all(~isnan(options.Shift)) && numel(options.Shift)~=numel(options.Nsegments)
        error('segmentize:shiftNsegments','The number of shifts is not equal to the number of segment numbers!');
    end
elseif any(isnan(options.Nsegments))
    if all(~isnan(options.Segsize)) && all(~isnan(options.Shift)) && numel(options.Shift)~=numel(options.Segsize)
        error('segmentize:shiftsegsize','The number of shifts is not equal to the number of segment lengths!');
    end
end

if any(~isnan(options.Segsize)) && ...
        (any(options.Segsize<0) || any(options.Segsize>size(X,options.Dim)))
    error('segmentize:segsize','The segment sizes cannot be negative or larger than the largest dimension!');
end

% Set the segment sizes and the shifts
if isnan(options.Segsize)
    if isnan(options.Shift)
        options.Segsize = calcnumbersegments(N,options.Nsegments,zeros(size(options.Nsegments)));
        options.Shift = cumprod([options.Segsize options.Nsegments(1:end-1)]);
    else
        options.Segsize = floor(N-sum(options.Shift.*options.Nsegments-options.Shift));
    end
    % Given the shifts, the segment sizes and the number of segments, determine the overlap
    overlap = calcoverlap([options.Segsize options.Nsegments(1:end-1)],options.Shift);
elseif isnan(options.Nsegments)
    if isnan(options.Shift), options.Shift = cumprod(options.Segsize); end
    % Given the shifts and the segment sizes, determine the overlap
    overlap = calcoverlap(options.Segsize,options.Shift);
    % Given the segment sizes and the overlap, determine the number of
    % segments
    options.Nsegments = calcnumbersegments(N,options.Segsize,overlap);
end

if any(isnan(options.Nsegments)) || any(isinf(options.Nsegments)), error('segmentize:Nsegmentsinf','The number of segments is not a natural number!'); end
if any(options.Nsegments<0), error('segmentize:Nsegmentsneg','The shift or segment size parameters have caused for negative segment sizes!'); end
if any(isnan(options.Segsize)) || any(isinf(options.Segsize)), error('segmentize:segsizeinf','The segment sizes are not a natural number!'); end
if any(options.Segsize<0), error('segmentize:Nsegmentsneg','The shift or segment size parameters have caused for a first segment size of %d, but this should not be negative!',options.Segsize); end

options.Segsize = [options.Segsize options.Nsegments(1:end-1)];
options.Nsegments = options.Nsegments(end);
size_segmentize = [options.Segsize options.Nsegments];

if options.UseAllSamples && N~=prod(size_segmentize)
    error('segmentize:notusingallsamples','The given values do not correspond with the number of samples!');
end

% Determine the storage needed
whosX = whos('X');
totalbytes = whosX.bytes/size(X,options.Dim)*prod(size_segmentize);
if ischar(options.Full) && strcmp(options.Full,'auto')
    if totalbytes > options.FullLimit*2^30
        options.Full = false;
        warning('segmentize:fullLimit',['The fullLimit has been reached, ',...
            'and instead of the dense tensor, the efficient ',...
            'representation is returned. The dense tensor can be obtained ',...
            'by setting the ''full'' option to true.']);
    else options.Full = true;
    end
end

repermorder = [options.Order+(1:options.PermToFirst-1) 1:options.Order options.Order+options.PermToFirst:numel([size_segmentize,size_other])];

if nargout>1 || ~options.Full
    % Construct the efficient representation
    sstruct = struct;
    sstruct.dim = options.Dim;
    sstruct.segsize = options.Segsize;
    sstruct.nsegments = options.Nsegments;
    sstruct.order = options.Order;
    sstruct.shift = options.Shift;
    sstruct.val = Xmat;
    sstruct.type = 'segment';
    sstruct.subsize.segment = size_segmentize;
    sstruct.subsize.other = size_other;
    sstruct.ispermuted = options.PermToFirst~=options.Dim;
    sstruct.repermorder = repermorder;
    sstruct.size = [size_segmentize size_other];
    sstruct.size = [size_segmentize size_other ...
        ones(1,numel(sstruct.repermorder)-numel([size_segmentize size_other]))];
    sstruct.size = sstruct.size(sstruct.repermorder);
end

reducehankel = false;
reducereshape = false;
if all(options.Shift == ones(1,options.Order-1))
    reducehankel = true;
else
    if ~isdefault.('Nsegments')
        if all(options.Shift == cumsum([options.Nsegments options.Segsize(end:-1:2)])-(1:numel(options.Segsize))+1)
            reducereshape = true;
        end
    else
        if all(options.Shift == cumsum(options.Segsize)-(1:numel(options.Segsize))+1)
            reducereshape = true;
        end
    end
end

if options.Full
    if reducehankel
        % Equivalent to Hankelization
        ind = cumsum(options.Segsize)-(1:numel(options.Segsize))+1;
        S = hankelize(Xmat,'Order',options.Order,'ind',ind);
        S = reshape(S,[size_segmentize size_other]);
    elseif reducereshape
        % Equivalent to basic reshape
        S = reshape(Xmat(1:prod(size_segmentize),:),[size_segmentize size_other]);
    else
        % Overlap present
        tmp = (1:size(Xmat,1));
        S = segvector(tmp,options.Segsize,options.Nsegments,overlap);
        S = reshape(Xmat(S,:),[size_segmentize size_other]);
    end
    
    % Permute the data to the specific dimensions
    if options.PermToFirst~=1
        S = permute(S,repermorder);
    end
else
    % Return the efficient representation
    S = sstruct;
end

end

function c = ccd(x,y,order)
% Closest common divisor: searches closest divisor to y of x, with y
% default sqrt(x)
if nargin<3, order = 2; end
if nargin<2, y = x.^(1/order); end
if isempty(y), y = x^(1/order); end
if abs((ceil(y)-y)/y) < 1e2*eps, tmp  = ceil(y);
else tmp = floor(y); end
if order==2, c = tmp;
else c = [tmp*ones(1,order-2) floor(x/tmp^(order-1))];
end
end

function d = calcnumbersegments(N,Segsizes,overlap)
% Determine number of segments in function of number of samples, segment
% sizes and overlap
if numel(Segsizes)~=numel(overlap), error('The number of segment sizes and the number of overlap parameters do not match!'); end
L = calcsamplesneeded(Segsizes,overlap);
d = floor((N-overlap(end))/(L-overlap(end)));
end

function S = segvector(v,Segsizes,Nsegmentsments,overlap)
% Determine the segmentation vectors
if numel(Segsizes)~=numel(overlap), error('The number of segment sizes and the number of overlap parameters do not match!'); end
if ~isvector(v), error('v should be a vector!'); end
v = v(:);
if numel(Segsizes)==1
    % Lowest level
    S = zeros([Segsizes Nsegmentsments]);
    for i = 1:Nsegmentsments
        idx = (i-1)*Segsizes-(i-1)*overlap+1:i*Segsizes-(i-1)*overlap;
        S(:,i) = v(idx,:);
    end
    S = reshape(S,[],1);
else
    S = zeros([prod(Segsizes) Nsegmentsments]);
    counter = 1;
    for i = 1:Nsegmentsments
        LL = calcsamplesneeded(Segsizes,overlap);
        idx = counter:counter+LL-1;
        tmp = v(idx,:);
        sv = segvector(tmp,Segsizes(1:end-1),Segsizes(end),overlap(1:end-1));
        S(:,i) = sv;
        counter = counter+LL-overlap(end);
    end
    if isempty(Nsegmentsments) || Nsegmentsments==0, S = [];
    else S = reshape(S,[],1);
    end
end
end

function l = calcsamplesneeded(Segsizes,overlap)
% Determine the samples needed, given the segment sizes and the overlap
NN = [1 cumprod(Segsizes(end:-1:1))];
l = NN(end);
for i = 1:numel(NN)-2
    l = l+overlap(i)*(-NN(end-i)+NN(end-i-1));
end
end

function overlap = calcoverlap(Segsize,shifts)
% Determine the overlap, given the segment sizes and the shifts
if numel(Segsize)~=numel(shifts)
    error('segmentize:myoverlap',['The number of segment sizes and the ',...
        'number of shifts do not match.']);
end
L = zeros(size(Segsize));
L(1) = Segsize(1);
for i = 2:numel(Segsize)
    L(i) = L(i-1)  + shifts(i-1)*(Segsize(i)-1);
end
overlap = L-shifts;
end