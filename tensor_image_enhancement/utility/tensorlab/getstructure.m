function structure = getstructure(varargin)
%GETSTRUCTURE determines the type of a tensor.
%   STRUCTURE = GETSTRUCTURE(T) detects the type or structure of the given
%   tensor T. STRUCTURE is a string containing the type. Currently, the
%   following structures are detected: full, incomplete, sparse, cpd, btd,
%   lmlra, tt, hankel and loewner.
%
%   STRUCTURE = GETSTRUCTURE(U,S) detects the structure of an LMLRA.
%
%   See also isvalidtensor.

% Author(s):   Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%              Otto Debals         (Otto.Debals@esat.kuleuven.be)
%              Lieven De Lathauwer (Lieven.DeLathauwer@esat.kuleuven.be)
%
% Version History:
% - 2015/12/16   NV      Initial version

if nargin == 1
    if isnumeric(varargin{1})
        structure = 'full';
    elseif iscell(varargin{1})
        if all(cellfun(@isnumeric, varargin{1}))
            if all(cellfun(@ndims, varargin{1})==2)
                structure  = 'cpd';
            else
                structure =  'tt';
            end
        elseif all(cellfun(@iscell, varargin{1}))
            structure = 'btd';
        elseif length(varargin{1}) > 1 && length(varargin{1})==2 && ...
                iscell(varargin{1}{1}) && isnumeric(varargin{1}{2}) 
            structure = 'lmlra';
        else
            error('getstructure:unknown', 'Unknown structure');
        end
    elseif isstruct(varargin{1})
        if isfield(varargin{1},'incomplete') && varargin{1}.incomplete
            structure = 'incomplete';
        elseif isfield(varargin{1},'sparse') && varargin{1}.sparse
            structure = 'sparse';
        elseif isfield(varargin{1},'type')
            structure = varargin{1}.type;
        else
            error('getstructure:unknown','Unknown structure');
        end
    else
        error('getstructure:unknown','Unknown structure');
    end
elseif nargin == 2
    if iscell(varargin{1}) && isnumeric(varargin{2})
        structure = 'lmlra';
    else
        error('getstructure:unknown', 'Unknown structure');
    end
else
    error('getstructure:unknown', 'Unknown structure');
end
end
