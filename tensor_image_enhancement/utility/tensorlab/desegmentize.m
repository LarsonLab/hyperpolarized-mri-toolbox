function X = desegmentize(S,varargin)
%DESEGMENTIZE Recover segmented signal(s).
%   X = DESEGMENTIZE(S) converts a matrix (or tensor) S of size
%   I1xI2x...xIN into a column vector X of length I1*I2*...*IN by unfolding
%   the matrix (tensor).
%
%   X = DESEGMENTIZE(S,'Order',K) converts a N-th order tensor of size
%   I1xI2x...xIN into a tensor of size (I1*I2*...*IK) x I_(K+1) x ... x
%   I_N by unfolding the first K modes.
%
%   X = DESEGMENTIZE(S,'Dims',n) converts an N-th order tensor by unfolding
%   the modes indicated by n. The resulting tensor has order N-length(n)+1.
%
%   X = DESEGMENTIZE(...,'Dim',d,'Order',K) is equivalent to
%   DESEGMENTIZE(S,'Dims',d:d+K-1).
%
%   DESEGMENTIZE(S,'key',value,...) or DESEGMENTIZE(S,options) can be used to
%   pass the following options:
%
%   - Shift:         Shift is a vector of length K-1, indicating the shift
%                    between consecutive segments. If S contains
%                    overlapping segments, this should be indicated with
%                    the Shift option.
%
%   - Method:        Indicates the desegmentization method:
%                    - If 'fibers', specific fibers (or part of fibers) are
%                      extracted.
%                    - If 'mean' or @mean (default), the overlapping
%                      parts are averaged.
%                    - If a function handle (such as @median), this
%                      function is applied on the overlapping parts.
%                    Alternatives of 'mean' can only be used when S is a
%                    full tensor.
%
%   - PermToDim:     Permutes the result such that the PermToDim is the detensorized
%                    mode. By default: Dims(1) or Dim.
%
%   - L:             When S is a polyadic representation, L is an array
%                    determining the number of columns of the factor
%                    matrices for each desegmentation. By default, L is
%                    equal to the total number of rank-1 terms.
%
%   See also segmentize, dehankelize, deloewnerize, dedecimate

%   Authors: Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
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
%   - 2015/11/18   OD      Optimized version

p = inputParser();
p.addOptional('Dim',1);
p.addOptional('Dims',NaN);
p.addOptional('Order',2);
p.addOptional('PermToDim',NaN);
p.addOptional('Shift',NaN);
p.addOptional('Rank',NaN);
p.addOptional('Method',@mean);
p.KeepUnmatched = false;
p.parse(varargin{:});
options = p.Results;
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
        if options.Dims+options.Order-1>getorder(S)
            error('desegmentize:dimsorder','The given order is not consistent with the dimensions of S!');
        end
        if numel(options.Dims)==1
            options.Dims = options.Dims:options.Dims+options.Order-1;
        end
    else
        % order is not set
        if any(options.Dims>getorder(S))
            error('desegmentize:dims','The given desegmentation dimensions are not consistent with the dimensions of S!');
        end
        options.Order = numel(options.Dims);
    end
    
    if numel(options.Dims)~=options.Order
        error('desegmentize:dimsandorder',...
            'The number of detensorized dimensions should be equal to the order!');
    end
else
    % dims is not set
    
    if ~isdefault.('Dim')
        % dim is set
        if options.Dim+options.Order-1>getorder(S)
            error('desegmentize:dim','The given dehankelization dimensions are not consistent with the dimensions of S!');
        end
    end
    
    if ~isdefault.('Order') && options.Order>getorder(S)
        error('desegmentize:order','The given order is not consistent with the dimensions of S!');
    end
    
    options.Dims = options.Dim:(options.Dim+options.Order-1);
end

if ~strcmp(getstructure(S),'cpd') && ~isdefault.('Rank')
    error('dehankelize:rank','The rank option is not supported when S is not a CPD!');
end

switch getstructure(S)
    case 'full'
        % Given a full tensor
        
        sh = size(S);
        size_segment = sh(options.Dims);
        segsize = size_segment(1:end-1);
        nsegs = size_segment(end);
        
        Smat = tens2mat(S,options.Dims);
        
        if isdefault.('Shift')
            options.Shift = cumsum(segsize)-(1:numel(segsize))+1;
        end
        
        overlap = cumsum(segsize)-options.Shift-(1:numel(options.Shift))+1;
        
        if strcmp(options.Method,'fibers')
            options.Method = @(x) x(1);
        elseif ischar(options.Method)
            options.Method = str2func(options.Method);
        end
            
        [sampleind,counts] = segvector(1:prod(size_segment),segsize,nsegs,overlap,zeros(prod(size_segment),1));
        N = sampleind(end);
        X = zeros(N,size(Smat,2));
        if isequal(options.Method,@mean), method = @sum;
        else method = options.Method;
        end
        
        for i = 1:size(Smat,2)
            X(:,i) = accumarray(sampleind,Smat(:,i),[],method);
        end
        
        if isequal(options.Method,@mean)
            X = bsxfun(@rdivide,X,counts(1:size(X,1)));
        end
        
        sx = sh; sx(options.Dims) = NaN;
        sx(options.Dims(1)) = N;
        sx(isnan(sx)) = [];
        signaldim = options.Dims(1)-sum(options.Dims(1)>options.Dims(2:end));
        X = mat2tens(X,sx,signaldim);
        
    case 'segment'
        % Given an efficient segmentized representation

        if ~isdefault.('Dims') || ~isdefault.('Order') || ~isdefault.('Shift')
            error('desegmentize:segment','The dims, order and shifts are used from the segment field')
        end
        if ~isvector(S.val)
            X = reshape(S.val,[size(S.val,1) S.subsize.other]);
            newrepermorder = S.repermorder;
            newrepermorder(newrepermorder > 1 & newrepermorder <= S.order) = [];
            newrepermorder(newrepermorder > S.order) = newrepermorder(newrepermorder>S.order)-S.order+1;
            X = permute(X,newrepermorder);
            
        else
            if S.dim==1, X = S.val;
            elseif S.dim==2, X = S.val.';
            else X = permute(S.val,[S.dim 1:S.dim]);
            end
        end
        if S.ispermuted, signaldim = 1;
        else signaldim = S.dim;
        end
    case 'cpd'
        % Given a representation of a tensor in rank-1 terms
        
        if any(cellfun('size', S, 2) ~= size(S{1},2))
            print('For a CPD size(S{n},2) should be R for all n\n');
            return;
        end
        if ~isdefault.('Rank')
            if any(options.Rank<1)
                error('desegmentize:rank1','The different ranks should be larger than 1!')
            end
            if sum(options.Rank)~=size(S{1},2)
                error('desegmentize:ranksum','The sum of the ranks should be equal the rank of the factor matrices!');
            end
        else
            options.Rank = size(S{1},2);
        end
        rank = options.Rank;
        ranke = cumsum([1 rank]);
        idx = find(strncmpi('Rank',varargin,3) | strncmpi('Rank',varargin,4));
        varargin(idx:idx+1) = [];
        if numel(rank)~=1
            Xt = cell(1,numel(rank));
            for r = 1:numel(rank)
                T = ful(cellfun(@(x) x(:,ranke(r):ranke(r+1)-1) ,S,'UniformOutput',false));
                Xt{r} = desegmentize(T,varargin{:});
            end
            X = cat(ndims(Xt{1}),Xt{:});
            X = permute(X,[1:options.Dims(1) ndims(X) options.Dims(1)+1:ndims(X)-1]);
        else 
            T = ful(S);
            X = desegmentize(T,varargin{:});
        end
    case {'incomplete','sparse','tt','btd','lmlragen'}
        T = ful(S);
        X = desegmentize(T,varargin{:});
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

end

function [S,counts] = segvector(v,segsizes,nsegments,overlap,counts)
% Calculate the (overlapping) index vectors
if nargin<5 && nargout==2, error('No counts in means no counts out!'); end
if numel(segsizes)~=numel(overlap), error('The number of segment sizes and the number of overlap parameters do not match!'); end
if ~isvector(v), error('v should be a vector!'); end
v = v(:);
if numel(segsizes)==1
    S = zeros([segsizes nsegments]);
    for i = 1:nsegments
        idx = (i-1)*segsizes-(i-1)*overlap+1:i*segsizes-(i-1)*overlap;
        S(:,i) = v(idx,:);
        if nargin>4, counts(idx) = counts(idx)+1; end
    end
    S = reshape(S,[],1);
else
    S = zeros([prod(segsizes) nsegments]);
    counter = 1;
    for i = 1:nsegments
        LL = leso(segsizes,overlap);
        idx = counter:counter+LL-1;
        tmp = v(idx,:);
        if nargin<5
            sv = segvector(tmp,segsizes(1:end-1),segsizes(end),overlap(1:end-1));
        else
            [sv,counttmp] = segvector(tmp,segsizes(1:end-1),segsizes(end),overlap(1:end-1),zeros(size(idx.')));
            counts(idx) = counts(idx) + counttmp;
        end
        S(:,i) = sv;
        counter = counter+LL-overlap(end);
    end
    if isempty(nsegments) || nsegments==0, S = [];
    else S = reshape(S,[],1);
    end
end
end

function l = leso(segsizes,overlap)
NN = [1 cumprod(segsizes(end:-1:1))];
l = NN(end);
for i = 1:numel(NN)-2
    l = l+overlap(i)*(-NN(end-i)+NN(end-i-1));
end
end
