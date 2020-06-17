function M = tens2mat(T,mode_row,mode_col)
%TENS2MAT Matricize a tensor.
%   M = tens2mat(T,mode_row,mode_col) matricizes a full or sparse tensor T into
%   a full or sparse matrix M of dimensions prod(size_tens(mode_row))-
%   by-prod(size_tens(mode_col)), where size_tens is equal to size(T). The
%   columns (rows) of M are obtained by fixing the indices of T corresponding to
%   mode_col (mode_row) and looping over the remaining indices in the order
%   mode_row (mode_col). E.g., if A and B are two matrices and T = cat(3,A,B),
%   then tens2mat(T,1:2,3) is the matrix [A(:) B(:)].
%
%   M = tens2mat(T,mode_row) matricizes a tensor T, where mode_col is
%   chosen as the sequence [1:ndims(T)]\mode_row.
%
%   M = tens2mat(T,[],mode_col) matricizes a tensor T, where mode_row is
%   chosen as the sequence [1:ndims(T)]\mode_col.
%
%   See also mat2tens.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/02   NV      Extension to sparse tensors
% - 2014/02/01   LS      Initial version
    
% Check arguments.
if isnumeric(T), size_tens = size(T); 
else size_tens = getsize(T); end
N = length(size_tens);
size_tens = [size_tens 1];
if nargin <= 2, mode_col = []; end
if isempty(mode_row) && isempty(mode_col)
    error('tens2mat:InvalidModes', ...
          'Either mode_row or mode_col must be non-empty.');
end
mode_row = mode_row(mode_row <= N); % >N is treated as singleton dimension.
mode_col = mode_col(mode_col <= N);
if isempty(mode_col), 
    mode_col = 1:N;
    mode_col(mode_row) = [];
end
if isempty(mode_row), 
    mode_row = 1:N;
    mode_row(mode_col) = [];
end
if isempty(mode_col), mode_col = N+1; end
if isempty(mode_row), mode_row = N+1; end

if isnumeric(T)
    % full tensor
    if any(mode_row(:).' ~= 1:length(mode_row)) || ...
            any(mode_col(:).' ~= length(mode_row)+(1:length(mode_col)))
        T = permute(T,[mode_row(:).' mode_col(:).']);
    end
    M = reshape(T,prod(size_tens(mode_row)),[]);

elseif isstruct(T) && isfield(T, 'sparse') && T.sparse
    % sparse tensor
    if ~isfield(T, 'sub'), T = fmt(T); end

    % create row indices
    if mode_row == N+1
        rows = 1;
        size_row = 1;
    elseif length(mode_row) == 1 && mode_row == 1 && isfield(T, 'matrix') && ...
            ~isempty(T.matrix)
        M = T.matrix;
        return;
    else 
        rows = double(T.sub{mode_row(1)});
        cumsize = cumprod([1 T.size(mode_row)]);
        for n = 2:length(mode_row)
            rows = rows + double(T.sub{mode_row(n)}-1)*cumsize(n);
        end
        size_row = cumsize(end);
    end
        
    % create column indices
    if mode_col == N+1
        cols = 1;
        size_col = 1;
    else 
        cols = double(T.sub{mode_col(1)});
        cumsize = cumprod([1 T.size(mode_col)]);
        for n = 2:length(mode_col)
            cols = cols + double(T.sub{mode_col(n)}-1)*cumsize(n);
        end
        size_col = cumsize(end);
    end 

    % sort for increased efficiency of sparse
    if length(rows) > 1
        [rows, j] = sort(rows);
        rows = double(rows);
        if length(cols) > 1, cols = cols(j); end
        vals = T.val(j);
    else 
        [cols, j] = sort(cols);
        vals = T.val(j);        
    end 
    M = sparse(rows(:), cols(:), vals(:), size_row, size_col);
else % other type
    error('tens2mat:invalidTensor', ['Only full and sparse tensors are supported. ' ...
                        'Use tens2mat(ful(T),mode_row,mode_col) for other types.']);
end
end
