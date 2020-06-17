function [D,dstruct] = decimate(X,varargin)
%DECIMATE Decimation of vectors, matrices or tensors.
%   D = DECIMATE(X) with a vector X of length N=IJ constructs a matrix of
%   size IxJ with the J consecutive segments of length I in its columns.
%
%   D = DECIMATE(X) with a matrix X constructs a matrix from the different
%   subsampled segments of every column, and stacks these matrices along
%   the third mode of a third-order tensor.
%
%   D = DECIMATE(X) with a tensor X constructs a matrix from the different
%   subsampled segments of every mode-1 fiber, and stacks them along the
%   higher-order modes of X.
%
%   [D,str] = DECIMATE(X,...) returns a structure str as an efficient
%   representation of the output tensor. str can be used for efficient
%   computations such as tensor decompositions and visualizations.
%
%   D = DECIMATE(...,'Order',K) constructs a tensor of order K from each
%   fiber by reshaping the data into a tensor of size I_1x...xI_K. The
%   values I_1,...,I_K are chosen such that they are as close together as
%   possible and such that as much data of X is used as possible. By default:
%   K=2.
%
%   D = DECIMATE(...,'Dim',n) constructs a matrix or higher-order tensor
%   (depending on the order) for every mode-n fiber. By default: n=1.
%
%   D = DECIMATE(X,...,'Full',full) avoids the construction of the full
%   tensor if Full is false. S is then the efficient representation
%   of the output tensor. This can be useful if the output tensor is
%   large. If Full is 'auto' (by default), the full tensor is returned if
%   the storage of the tensor would not exceed the value of the FullLimit
%   option; otherwise, the efficient representation is returned. If Full is
%   true, the full tensor is always returned.
%
%   DECIMATE(X,'key',value) or DECIMATE(X,options) can be used to
%   pass the following options:
%
%   - Nsamples:         Instead of the default values I_1,...,I_K-1, the
%                       values from Nsamples are used. I_K is chosen such
%                       that as much data of X is used as possible.
%                       Nsamples is a vector of length K-1. If the
%                       Subsample option is used (see further), an error is
%                       thrown.
%
%   - Subsample:        Instead of the default values I_2,...,I_K, the
%                       values from Subsample are used. I_1 is chosen such
%                       that as much data of X is used as possible.
%                       Subsample is a vector of length K-1. If the
%                       Nsamples option is used as well, an error is thrown.
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
%   See also dedecimate, hankelize, loewnerize, segmentize

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
p.addOptional('Subsample', NaN);
p.addOptional('Nsamples',NaN);
p.addOptional('Order', 2);
p.addOptional('Dim', 1);
p.addOptional('PermToFirst',false);
p.addOptional('Full','auto');
p.addOptional('Shift',NaN);
p.addOptional('UseAllSamples',false);
p.addOptional('FullLimit',1);
p.KeepUnmatched = false;
p.parse(varargin{:});

% Processing inputs
isdefault = cell2struct(cellfun(@(x) any(strcmp(x,p.UsingDefaults)),p.Parameters,...
    'UniformOutput',false),p.Parameters,2);

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

if ~isdefault.('Subsample') && ~isdefault.('Nsamples')
    error('Only one of the options subsample and nsamples is allowed to be used at the same time')
end
if ~isdefault.('Subsample')
    % Revert the subsample option values to use segmentize
    options.Segsize = options.Subsample(end:-1:1);
end
options = rmfield(options,'Subsample');

if ~isdefault.('Nsamples')
    % Revert the nsamples option values to use segmentize
    options.Nsegments = options.Nsamples(end:-1:1);
end
options = rmfield(options,'Nsamples');

for i = 1:numel(p.UsingDefaults)
    if ~any(strcmp(p.UsingDefaults{i},{'Subsample','Nsamples'}))
        options = rmfield(options,p.UsingDefaults{i});
    end
end
[D,dstruct] = segmentize(X,options);

% Correct the efficient representation
dstruct.subsample = dstruct.segsize(end:-1:1);
dstruct = rmfield(dstruct,'segsize');
dstruct.nsamples = dstruct.nsegments;
dstruct = rmfield(dstruct,'nsegments');
dstruct.type = 'decimate';
dstruct.subsize.segment = dstruct.subsize.segment(end:-1:1);

[~,idx] = sort(dstruct.repermorder);
idx = idx(1:dstruct.order);
dstruct.size(idx) = dstruct.size(idx(end:-1:1));

if strcmp(getstructure(D),'segment')
    D = dstruct;
else
    v = 1:numel(dstruct.size);
    v(idx) = v(idx(end:-1:1));
    if any(v~=1:numel(v))
        D = ipermute(D,v);
    end
end