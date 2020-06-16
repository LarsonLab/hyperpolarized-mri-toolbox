function c = dcov(x,varargin)
%DCOV Covariance matrices along specific dimensions.
%   C = DCOV(X) with a vector X with N observations computes the variance of X.
%
%   C = DCOV(X) with a matrix X with N rows (the observations) and R
%   columns (the variables) constructs the covariance matrix of X of size
%   RxR.
%
%   C = DCOV(X) with a tensor X constructs a covariance matrix for every
%   slice of X assuming that the first dimension represents the
%   observations and the second dimension represents the variables. The
%   covariance matrices are stacked along the other dimensions.
%
%   C = DCOV(...,'Dims',d) calculates the covariance matrices by assuming
%   that d(1) represents the dimension with the observations and that d(2)
%   represents the dimension with the variables.
%
%   DCOV(X,'key',value) or DCOV(X,options) can be used to pass the
%   following options:
%   
%   - Unbiased: If true (default), the covariance matrices are normalized
%               by N-1. Otherwise, a normalization with N is used.
%
%   - PermToDims: The covariance matrices or covariance tensors are inserted
%               along the PermToDims modes, instead of at the original
%               tensorized dimensions. Default: Dims.
%
%   - Missing:  If 'includenan', the output contains NaN if the input
%               contains NaN or missing elements. If 'omitrows', all of the
%               rows of X containing NaN values are omitted. If
%               'partialrows', each element C(I,J) is computed separately,
%               based only on the columns I and J of X. Rows are only
%               omitted if they contain NaN values in column I or J of X.
%
%   See also cov, hankelize, loewnerize, segmentize, decimate

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
%   - 2016/02/13    OD      Initial version


p = inputParser();
p.addOptional('Dims', [1 2]);
p.addOptional('Unbiased', true);
p.addOptional('PermToDims',NaN);
p.addOptional('Missing','includenan');
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;

% Processing inputs
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);

nondims = 1:ndims(x); nondims(options.Dims) = [];
Xp = permute(x,[options.Dims nondims]);
s = size(Xp);

if numel(s)>2, size_other = s(3:end);
else size_other = [];
end

Xp = reshape(Xp,[s(1) s(2) prod(size_other)]);

c = zeros([s(2) s(2) prod(size_other)]);
for i = 1:size(Xp,3)
    if strcmp('missing',p.UsingDefaults)
        c(:,:,i) = cov(Xp(:,:,i),~options.Unbiased,options.Missing);
    else
        c(:,:,i) = cov(Xp(:,:,i),~options.Unbiased);
    end
end
c = reshape(c,[s(2) s(2) size_other]);

if isdefault.('PermToDims'), options.PermToDims = options.Dims; end
if islogical(options.PermToDims)
    if options.PermToDims, options.PermToDims = [1 2];
    else options.PermToDims = options.Dims;
    end
end
if numel(options.PermToDims)==1
    options.PermToDims = options.PermToDims+(0:1);
end

if any(options.PermToDims~=[1 2])
    [mindim,mindimidx] = min(options.PermToDims);
    [maxdim,maxdimidx] = max(options.PermToDims);
    v = [3:mindim+1 mindimidx mindim+2:maxdim maxdimidx maxdim+1:ndims(x)];
    c = permute(c,v);
end