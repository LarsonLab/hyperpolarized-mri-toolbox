function sz = getsize(T, dim)
%GETSIZE Dimensions of a tensor.
%   SZ = GETSIZE(T) returns the dimensions of a tensor T. The tensor can be
%   given as an array of numerical values, a struct, or a decomposition (a
%   cell). In the last case, the type of the decomposition is inferred from the
%   structure of the cell.
%
%   SZ = GETSIZE(T,DIM) returns the dimensions of tensor T, specified by DIM.
%   Multiple dimensions can be requested at once. 
%
%   See also getorder.

% Authors: Nico Vervliet        (Nico.Vervliet@esat.kuleuven.be)
%          Otto Debals          (Otto.Debals@esat.kuleuven.be)
%          Lieven De Lathauwer  (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2015/12/06   NV      Initial version

try
    switch getstructure(T)
        case 'full'
            sz = size(T);
        case {'incomplete', 'sparse'}
            if isfield(T, 'size')
                sz  = T.size;
            else
                error('getsize:T', 'T is has an invalid structure');
            end
        case 'cpd'
            sz = cellfun('size', T, 1);
        case 'lmlra'
          sz = cellfun('size', T{1}, 1);
        case 'btd'
            sz = cellfun('size', T{1}(1:end-1), 1);
        case 'tt'
            sz = cellfun('size',T,2);
            sz(1) = size(T{1},1);
        case 'hankel'
            sz = T.size;
        case 'loewner'
            sz = T.size;
        case 'segment'
            sz = T.size;
        case 'decimate'
            sz = T.size;
        otherwise
            error('getsize:invalidTensor', 'T has an unknown structure');
    end
    if nargin >= 2
        if any(dim <= 0) || any(dim > length(sz))
            error('getsize:dim', 'dim should be > 0 and <= getorder(T)');
        end
        sz = sz(dim);
    end
catch e
    if strcmpi(e.identifier, 'getstructure:unknown')
        error('getsize:invalidTensor', 'T has an unknown structure');
    else
        rethrow(e);
    end
end
