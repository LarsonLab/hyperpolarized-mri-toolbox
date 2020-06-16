function X = dedecimate(D,varargin)
%DEDECIMATE Recover decimated signal(s).
%   X = DEDECIMATE(D) converts a matrix (or tensor) D of size I1xI2x...x IN
%   into a column vector X of length I1*I2*...*IN by unfolding the matrix
%   (tensor)
%
%   X = DEDECIMATE(D,'Order',K) converts a N-th order tensor of size
%   I1xI2x...xIN into a tensor of size (I1*I2*...*IK)xI_(K+1)x...xI_N by
%   unfolding the first K modes.
%
%   X = DEDECIMATE(D,'Dims',n) converts an N-th order tensor by unfolding
%   the modes indicated by n. The resulting tensor has order N-length(n)+1.
%
%   X = DEDECIMATE(...,'Dim',d,'Order',K) is equivalent to
%   DEDECIMATE(D,'Dims',d:d+K-1).
%
%   DEDECIMATE(D,'key',value,...) or DEDECIMATE(D,options) can be used to
%   pass the following options:
%
%   - Shift:         Shift is a vector of length K-1, indicating the shift
%                    between consecutive subsampled segments. If D contains
%                    overlapping segments, this should be indicated with
%                    the Shift option.
%
%   - Method:        Indicates the dedecimation method:
%                    - If 'fibers', specific fibers (or part of fibers) are
%                      extracted.
%                    - If 'mean' or @mean (default), the overlapping
%                      parts are averaged.
%                    - If a function handle (such as @median), this
%                      function is applied on the overlapping parts.
%                    Alternatives of 'mean' can only be used when D is a
%                    full tensor.
%
%   - PermToDim:     Permutes the result such that the PermToDim is the detensorized
%                    mode. By default: Dims(1) or Dim.
%
%   - Rank:          When D is a polyadic representation, L is an array
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
        error('dedecimate:dimanddims','Using both the ''Dim'' and ''Dims'' option arguments is invalid!');
    end
    
    if ~isdefault.('Order')
        % order is set
        if options.Dims+options.Order-1>getorder(D)
            error('dedecimate:dimsorder','The given order is not consistent with the dimensions of S!');
        end
        if numel(options.Dims)==1
            options.Dims = options.Dims:options.Dims+options.Order-1;
        end
    else
        % order is not set
        if any(options.Dims>getorder(D))
            error('dedecimate:dims','The given dedecimation dimensions are not consistent with the dimensions of S!');
        end
        options.Order = numel(options.Dims);
    end
    
    if numel(options.Dims)~=options.Order
        error('dedecimate:dimsandorder',...
            'The number of detensorized dimensions should be equal to the order!');
    end
else
    % dims is not set
    
    if ~isdefault.('Dim')
        % dim is set
        if options.Dim+options.Order-1>getorder(D)
            error('dedecimate:dim','The given dedecimated dimensions are not consistent with the dimensions of !');
        end
    end
    
    if ~isdefault.('Order') && options.Order>getorder(D)
        error('dedecimate:order','The given order is not consistent with the dimensions of S!');
    end
    
    options.Dims = options.Dim:(options.Dim+options.Order-1);
end

if ~strcmp(getstructure(D),'cpd') && ~isdefault.('Rank')
    error('dedecimate:rank','The rank option is not supported when S is not a CPD!');
end

switch getstructure(D)
    case 'full'
        % Given a full tensor
        
        decpermorder = 1:ndims(D);
        decpermorder(options.Dims) = options.Dims(end:-1:1);
        D = permute(D,decpermorder);
        desegmentizeoptions = struct;
        desegmentizeoptions.Dims = options.Dims;
        if ~isdefault.('Shift'), desegmentizeoptions.Shift = options.Shift; end
        if ~isdefault.('PermToDim'), desegmentizeoptions.PermToDim = options.PermToDim; end
        X = desegmentize(D,desegmentizeoptions);
    case 'decimate'
        % Given an efficient decimated representation
        
        if ~isdefault.('Dims') || ~isdefault.('Order') || ~isdefault.('Shift')
            error('dedecimate:decimate','The dims, order and shifts are used from the decimate field')
        end
        
        D.segsize = D.subsize.segment(end:-1:2);
        D = rmfield(D,'subsample');
        D.type = 'segment';
        desegmentizeoptions = struct;
        if ~isdefault.('PermToDim'), desegmentizeoptions.PermToDim = options.PermToDim; end
        X = desegmentize(D,desegmentizeoptions);
        
    case 'cpd'
        % Given a representation of a tensor in rank-1 terms
        
        if any(cellfun('size', D, 2) ~= size(D{1},2))
            print('For a CPD size(S{n},2) should be R for all n\n');
            return;
        end
        if ~isdefault.('Rank')
            if any(options.Rank<1)
                error('dedecimate:rank1','The different ranks should be larger than 1!')
            end
            if sum(options.Rank)~=size(D{1},2)
                error('dedecimate:ranksum','The sum of the ranks should be equal the rank of the factor matrices!');
            end
        else
            options.Rank = size(D{1},2);
        end
        rank = options.Rank;
        ranke = cumsum([1 rank]);
        idx = find(strncmpi('Rank',varargin,3) | strncmpi('Rank',varargin,4));
        varargin(idx:idx+1) = [];
        if numel(rank)~=1
            Xt = cell(1,numel(rank));
            for r = 1:numel(rank)
                T = ful(cellfun(@(x) x(:,ranke(r):ranke(r+1)-1) ,D,'UniformOutput',false));
                Xt{r} = desegmentize(T,varargin{:});
            end
            X = cat(ndims(Xt{1}),Xt{:});
            X = permute(X,[1:options.Dims(1) ndims(X) options.Dims(1)+1:ndims(X)-1]);
        else 
            T = ful(D);
            X = dedecimate(T,varargin{:});
        end
    case {'incomplete','sparse','tt','btd','lmlragen'}
        T = ful(D);
        X = dedecimate(T,varargin{:});
    otherwise
        error('Structure not supported!')
end
end
