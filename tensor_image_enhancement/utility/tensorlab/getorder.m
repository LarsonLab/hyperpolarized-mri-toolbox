function N = getorder(T, ~)
%GETORDER Order of a tensor.
%   N = GETORDER(T) determines the order of a tensor T. T can be a numerical
%   array of values, a struct or a decomposition (a cell). In the last case,
%   the type of the decomposition is inferred from the structure of the cell.
%
%   N = GETORDER(U, S) determines the order for a LMLRA in which U are the
%   factor matrices U and S the core tensor.
%
%   See also getsize.

% Authors: Nico Vervliet        (Nico.Vervliet@esat.kuleuven.be)
%          Otto Debals          (Otto.Debals@esat.kuleuven.be)
%          Lieven De Lathauwer  (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2015/12/06   NV      Initial version
% - 2015/12/10   OD      Added Hankel/Loewner

try
    switch getstructure(T)
        case 'full'
            N = ndims(T);
        case {'incomplete', 'sparse'}
            if isfield(T, 'size')
                N  = length(T.size);
            else
                error('getorder:T', 'T is has an invalid structure');
            end
        case 'cpd'
            N = length(T);
        case 'lmlra'
            if nargin < 2
                N = length(T{1});
            else
                N = length(T);
            end
        case 'btd'
            N = length(T{1})-1;
        case 'tt'
            N = length(T);
        case 'hankel'
            if isfield(T, 'size')
                N = length(T.size);
            else
                error('getorder:T', 'T is has an invalid structure');
            end
        case 'loewner'
            if isfield(T, 'size')
                N = length(T.size);
            else
                error('getorder:T', 'T is has an invalid structure');
            end
        case 'segment'
            if isfield(T, 'size')
                N = length(T.size);
            else
                error('getorder:T', 'T is has an invalid structure');
            end
        case 'decimate'
            if isfield(T, 'size')
                N = length(T.size);
            else
                error('getorder:T', 'T is has an invalid structure');
            end
        otherwise
            error('getorder:invalidTensor', 'T has an unknown structure');
    end
catch e
    if strcmpi(e.identifier, 'getstructure:unknown')
        error('getorder:invalidTensor', 'T has an unknown structure');
    else
        rethrow(e);
    end
end
