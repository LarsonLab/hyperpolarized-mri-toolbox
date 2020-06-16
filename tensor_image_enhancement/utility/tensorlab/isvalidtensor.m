function isvalid = isvalidtensor(T, varargin)
%ISVALIDTENSOR Check if the representation of a tensor is correct.
%   ISVALIDTENSOR(T) returns true if the tensor is valid. A tensor is valid
%   if its type can be determined using GETSTRUCTURE(T) and if it is
%   formatted correctly. Some examples:
%      - CPD: all entries in the cell should be matrices with R columns.   
%      - Incomplete: the indices should be within the bounds determined by
%        the size of the tensor, the number of values should be equal to the
%        number of indices, no indices can be defined twice. etc. 
%   
%   ISVALIDTENSOR(T, printerrors) prints a message to state why a tensor is
%   invalid (if it is invalid) if printerrors = true.
%
%   ISVALIDTENSOR(T, printerrors, type) does not check the type of the tensor
%   using GETSTRUCTURE(T), but uses the given type. 
%
%   ISVALIDTENSOR(U,S), ISVALIDTENSOR(U,S,printerrors) and
%   ISVALID(U,S,printerrors,type) can be used for the LMLRA format. 
%   

% Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)     
%
% Version History:
% - 2015/12/21   NV      Initial version

% Parse input options
try 
    if isempty(varargin)
        printerrors = false;
        type = getstructure(T);
    elseif length(varargin) == 1
        if islogical(varargin{1})
            printerrors = varargin{1};
        else
            T = {T,varargin{1}};
            printerrors = false;
        end
        type = getstructure(T);
    elseif length(varargin) == 2
        if islogical(varargin{1})
            printerrors = varargin{1};
            type = varargin{2};
        else
            T = {T,varargin{1}};
            printerrors = varargin{2};
            type = getstructure(T);
        end
    elseif length(varargin) == 3
        T = {T, varargin{1}};
        printerrors = varargin{2};
        type = varargin{3};
    end
catch e 
    if strcmpi(e.identifier, 'getstructure:unknown')
        error('isvalidtensor:unknown', ['Not enough information available to ' ...
                            'determine the type of T.']);
    else
        rethrow(e);
    end
end

isvalid = false;
switch type
  case 'full'
    isvalid = true;
    return;
  case {'incomplete', 'sparse'}
    if ~isfield(T, 'size') || (~isfield(T, 'ind') && ~isfield(T, 'sub')) ...
            || ~isfield(T, 'val')
        print(['T should have at least the fields size, ind or sub, ' ...
               'and val.\n']);
        return;
    else
        if isfield(T, 'ind')
            if any(diff(sort(T.ind))==0)
                print('T has duplicate indices (in T.ind).\n');
                return;
            end
            if any(T.ind <= 0) || any(T.ind > prod(T.size))
                print(['T.ind has out-of-range indices (<= 0 or > ' ...
                       'prod(T.size))\n']);
                return;
            end
            len = length(T.ind);
        end
        if isfield(T, 'sub')
            try
                tmp = sub2ind(T.size, T.sub{:});
            catch
                print(['T.sub has out-of-range indices (<= 0 or > ' ...
                       'T.size(n) for some n)\n']);
                return;
            end
            if any(diff(sort(tmp))==0)
                print('T has duplicate indices (in T.sub).\n');
                return;
            end
            len = length(T.sub{1});
        end
        if isfield(T, 'ind') && isfield(T, 'sub')
            if length(T.ind) ~= length(tmp)
                print(['T.ind and T.sub{n} have not the same length for ' ...
                       'some n\n']);
                return;
            end
            if any(T.ind ~= tmp)
                print(['The indices for T.ind and T.sub do not ' ...
                       'correspond\n']);
                return;
            end
        end
        if length(T.val) ~= len
            print('T.val has not the same length as T.ind or T.sub\n');
            return;
        end
        if strcmpi(type, 'incomplete')
            if isfield(T, 'sparse') && T.sparse
                print(['T cannot be incomplete and sparse at the same ' ...
                       'time\n']);
                return
            end
        else
            if isfield(T, 'incomplete') && T.incomplete
                print(['T cannot be incomplete and sparse at the same ' ...
                       'time\n']);
                return
            end
        end
    end
  case 'cpd'
    if any(cellfun('size', T, 2) ~= size(T{1},2))
        print('For a CPD size(T{n},2) should be R for all n\n');
        return;
    end
  case 'lmlra'
    T{1} = T{1}(:)';
    if any(~cellfun(@ismatrix, T{1}))
        print('For an LMLRA T{1}{n} should be matrices for all n\n');
        return;
    end
    if find(cellfun('size', T{1}, 2) > 1, 1, 'last') > ndims(T{2}) || ...
            length(T{1}) < ndims(T{2})
        print('For an LMLRA length(T{1}) should be equal to ndims(T{2})\n');
        return;
    end
    if any(cellfun('size', T{1}(1:ndims(T{2})), 2) ~= size(T{2}))
        print(['For an LMLRA size(T{1}{n},2) should be size(T{2},n) for all ' ...
               'n\n']);
        return;
    end
  case 'btd'
    T = T(:)';
    % test per term correctness
    for r = 1:length(T)
        T{r} = T{r}(:)';
        if any(~cellfun(@ismatrix, T{r}(1:end-1)))
            print('For a BTD T{r}{n} should be matrices for all n=1:length(T{r})-1\n');
            return;
        end
        if find(cellfun('size', T{r}(1:end-1), 2) > 1, 1, 'last') > ...
                ndims(T{r}{end}) || length(T{r}) - 1 < ndims(T{r}{end})
            print(['For a BTD length(T{r}(1:end-1)) should be equal to ' ...
                   'ndims(T{r}{end}) for all r\n']);
            return;
        end
        if any(cellfun('size', T{r}(1:ndims(T{r}{end})), 2) ~= size(T{r}{end}))
            print(['For a BTD size(T{r}{n},2) should be size(T{r}{end},n) for n ' ...
                   '= 1:length(T{r})-1\n']);
            return;
        end
    end
    % test cross term correctness
    sizes = cell(1, length(T));
    for r = 1:length(T)
        sizes{r} = cellfun('size', T{r}(1:end-1), 1);
        sizes{r} = sizes{r}(1:find(sizes{r}>1,1,'last'));
    end
    if any(diff(cellfun(@length, sizes))~=0) || any(any(diff(cat(1,sizes{:}),1,1)~=0))
        print(['For a BTD size(T{r}{n},1) == size(T{s}{n},1) should hold ' ...
               'for all r,s=1:length(T), and n<=ndims(T{t}{end}) for all ' ...
               't\n']);
        return;
    end
  case 'tt'
    if ~ismatrix(T{1}) || ~ismatrix(T{end})
        print('For a TT, T{1} and T{end} should be matrices\n');
        return;
    end
    if any(cellfun(@ndims, T(2:end-1))>3)
        print(['For a TT, T{n} should be third order tensors for n = ' ...
               '1:length(T)-1\n']);
        return;
    end
    if size(T{1},2) ~= size(T{2},1)
        print('For a TT, size(T{1},2) should equal size(T{2},1)\n');
        return;
    end
    if any(cellfun('size', T(2:end-1), 3)~=cellfun('size', T(3:end), 1))
        print(['For a TT, size(T{n},3) should equal size(T{n+1},1) for n ' ...
               '= 2:length(T)-1\n']);
        return;
    end
  case 'hankel'
    % fields missing
    if ~isfield(T, 'size') || ~isfield(T, 'subsize') || ...
            ~isfield(T,'repermorder') || ~isfield(T,'ispermuted') || ...
            ~isfield(T,'val') || ~isfield(T,'dim') || ...
            ~isfield(T,'order') || ~isfield(T,'ind')
        print(['T should have at least the fields val, dim, order, ',...
               'ind, size, subsize, repermorder and ispermuted.\n']);
        return;
    end
    % no other/hankel subsize field
    if ~isfield(T.subsize,'hankel') || ~isfield(T.subsize,'other')
        print(['T.subsize should have at least the fields hankel and ',...
               'other.\n']);
        return;
    end
    % tensor for H.val
    if ~ismatrix(T.val)
        print('T.val should be a vector or a matrix.\n\')
        return;
    end
    % order < 2
    if T.order < 2
        print('T.order must be larger than or equal to 2.\n');
        return;
    end
    % number of ind != order-1
    if numel(T.ind)~=T.order-1
        print('The number of T.ind is not equal to T.order-1.\n');
        return;
    end
    % ind is not increasing
    if any(T.ind(2:end)<T.ind(1:end-1))
        print('T.ind should be increasing.\n');
        return;
    end
    % any ind is larger than the number of values
    if T.ind(end)>size(T.val,1)
        print('T.ind should be smaller than N.\n');
        return;
    end
    % number of size is smaller than order, or than order + 1
    if size(T.val,2)==1
        if numel(T.size) ~= T.order
            print(['If T.val is a vector, the dimension of T should '...
                   'equal T.order.\n']);
            return;
        end
    else
        if numel(T.size) <= T.order
            print(['If T.val is a matrix, the dimension of T should '...
                   'exceed T.order.\n']);
            return;
        end
    end
    % number of subsize.hankel elements is not equal to order
    if numel(T.subsize.hankel)~=T.order
        print('The number of elements in T.subsize.hankel should equal T.order.\n');
        return;
    end
    % prod(T.subsize.other) should equal size(T.val,2)
    if ~isempty(T.subsize.other) && prod(T.subsize.other)~=size(T.val,2)
        print('prod(T.subsize.other) is not equal to size(T.val,2).\n'); 
        return;
    end
    % sum(T.subsize.hankel)-order+1 should equal N
    if sum(T.subsize.hankel)-T.order+1 ~= size(T.val,1)
        print('sum(T.subsize.hankel)-T.order+1 should equal N.\n');
        return; 
    end
    % concatenation of subsize.hankel and subsize.other, and then
    % permute, is equal to size
    if T.ispermuted
        if any([T.subsize.hankel T.subsize.other]~=T.size)
            print(['T.subsize.other and T.subsize.hankel do not '...
                   'correspond with T.size.\n']);
            return;
        end
    else
        tmp = [T.subsize.hankel T.subsize.other];
        tmp = tmp(T.repermorder);
        if any(tmp~=T.size)
            print(['T.subsize.other and T.subsize.hankel do not '...
                   'correspond with T.size.\n']);
            return;
        end
    end
    % size and repermorder do not agree
    if T.ispermuted
        if any(T.repermorder~=1:numel(T.repermorder))
            print('If T.ispermuted is true, T.repermorder should equal 1:K.\n')
            return;
        end
    end
  case 'loewner'
    % fields missing
    if ~isfield(T, 'size') || ~isfield(T, 'subsize') || ...
            ~isfield(T,'repermorder') || ~isfield(T,'ispermuted') || ...
            ~isfield(T,'val') || ~isfield(T,'dim') || ...
            ~isfield(T,'order') || ~isfield(T,'ind') || ...
            ~isfield(T,'t') || ~isfield(T,'isequidistant')
        print(['T should have at least the fields val, dim, order, ',...
               'ind, size, subsize, repermorder, ispermuted, t and isequidistant.\n']);
        return;
    end
    % no other/loewner subsize field
    if ~isfield(T.subsize,'loewner') || ~isfield(T.subsize,'other')
        print(['T.subsize should have at least the fields loewner and ',...
               'other.\n']);
        return;
    end
    % tensor for H.val
    if ~ismatrix(T.val)
        print('T.val should be a vector or a matrix.\n\')
        return;
    end
    % order < 2
    if T.order < 2
        print('T.order must be larger than or equal to 2.\n');
        return;
    end
    % number of ind != order
    if numel(T.ind)~=T.order
        print('The number of T.ind is not equal to T.order.\n');
        return;
    end
    % number of indices is not equal to N
    if numel(cell2mat(T.ind))~=size(T.val,1)
        print('The number of indices should equal size(T.val,1).\n');
        return;
    end
    % indices are appearing multiple times
    u = unique(cell2mat(T.ind));
    if numel(u)~=size(T.val,1)
        print('Each index should only appear once across the index sets.\n');
        return; 
    end
    % index is larger than sample
    if u(end)>size(T.val,1)
        print('There are indices larger than size(T.val,1).\n');
        return;
    end
    if u(1)<1
        print('The indices should not be zero or negative.\n');
        return;
    end
    % number of size is smaller than order, or than order + 1
    if size(T.val,2)==1
        if numel(T.size) ~= T.order
            print(['If T.val is a vector, the dimension of T should '...
                   'equal T.order.\n']);
            return;
        end
    else
        if numel(T.size) <= T.order
            print(['If T.val is a matrix, the dimension of T should '...
                   'exceed T.order.\n']);
            return;
        end
    end
    % number of subsize.loewner elements is not equal to order
    if numel(T.subsize.loewner)~=T.order
        print('The number of elements in T.subsize.loewner should equal T.order.\n');
        return;
    end
    % prod(T.subsize.other) should equal size(T.val,2)
    if ~isempty(T.subsize.other) && prod(T.subsize.other)~=size(T.val,2)
        print('prod(T.subsize.other) is not equal to size(T.val,2).\n'); 
        return;
    end
    % sum(T.subsize.loewner) should equal N
    if sum(T.subsize.loewner)~= size(T.val,1)
        print('sum(T.subsize.loewner) should equal N.\n');
        return; 
    end
    % concatenation of subsize.loewner and subsize.other, and then
    % permute, is equal to size
    if T.ispermuted
        if any([T.subsize.loewner T.subsize.other]~=T.size)
            print(['T.subsize.other and T.subsize.loewner do not '...
                   'correspond with T.size.\n']);
            return;
        end
    else
        tmp = [T.subsize.loewner T.subsize.other];
        tmp = tmp(T.repermorder);
        if any(tmp~=T.size)
            print(['T.subsize.other and T.subsize.loewner do not '...
                   'correspond with T.size.\n']);
            return;
        end
    end
    % size and repermorder do not agree
    if T.ispermuted
        if any(T.repermorder~=1:numel(T.repermorder))
            print('If T.ispermuted is true, T.repermorder should equal 1:K.\n')
            return;
        end
    end
    % t is not a vector
    if ~isvector(T.t) || size(T.t,2)~=1
        print('t should be a column vector.\n');
        return;
    end
    % N and t do not agree
    if size(T.val,1) ~= numel(T.t)
        print('t should have size(T.val,1) elements.\n');
        return;
    end
  case 'segment'
    % fields missing
    if ~isfield(T, 'size') || ~isfield(T, 'subsize') || ...
            ~isfield(T,'repermorder') || ~isfield(T,'ispermuted') || ...
            ~isfield(T,'val') || ~isfield(T,'dim') || ...
            ~isfield(T,'order') || ~isfield(T,'segsize') || ...
            ~isfield(T,'nsegments') || ~isfield(T,'shift')
        print(['T should have at least the fields val, dim, order, ',...
               'ind, size, subsize, repermorder, ispermuted, segsize, nsegments and shift\n']);
        return;
    end
    % no other/segment subsize field
    if ~isfield(T.subsize,'segment') || ~isfield(T.subsize,'other')
        print(['T.subsize should have at least the fields segment and ',...
               'other.\n']);
        return;
    end
    % tensor for H.val
    if ~ismatrix(T.val)
        print('T.val should be a vector or a matrix.\n\')
        return;
    end
    % order < 2
    if T.order < 2
        print('T.order must be larger than or equal to 2.\n');
        return;
    end
    % number of segment sizes != order-1
    if numel(T.segsize)~=T.order-1
        print('The number of T.ind is not equal to T.order-1.\n');
        return;
    end
    % number of segment numbers != 1
    if numel(T.nsegments)~=1
        print('The number of segment numbers is not equal to 1.\n');
        return;
    end
    % number of shifts != T.order-1
    if numel(T.shift)~=T.order-1
        print('The number of shifts is not equal to T.order-1.\n'); 
        return;
    end
    % negative shifts
    if any(T.shift<=0)
        print('The shifts cannot be zero or negative.\n');
        return;
    end
    % number of size is smaller than order, or than order + 1
    if size(T.val,2)==1
        if numel(T.size) ~= T.order
            print(['If T.val is a vector, the dimension of T should '...
                   'equal T.order.\n']);
            return;
        end
    else
        if numel(T.size) <= T.order
            print(['If T.val is a matrix, the dimension of T should '...
                   'exceed T.order.\n']);
            return;
        end
    end
    % number of subsize.segment elements is not equal to order
    if numel(T.subsize.segment)~=T.order
        print('The number of elements in T.subsize.segment should equal T.order.\n');
        return;
    end
    % prod(T.subsize.other) should equal size(T.val,2)
    if ~isempty(T.subsize.other) && prod(T.subsize.other)~=size(T.val,2)
        print('prod(T.subsize.other) is not equal to size(T.val,2).\n'); 
        return;
    end
    % sum(T.subsize.segment)-order+1 should equal N
    if sum(T.subsize.segment)-T.order+1 ~= size(T.val,1)
        print('sum(T.subsize.segment)-T.order+1 should equal N.\n');
        return; 
    end
    % concatenation of subsize.segment and subsize.other, and then
    % permute, is equal to size
    if T.ispermuted
        if any([T.subsize.segment T.subsize.other]~=T.size)
            print(['T.subsize.other and T.subsize.segment do not '...
                   'correspond with T.size.\n']);
            return;
        end
    else
        tmp = [T.subsize.segment T.subsize.other];
        tmp = tmp(T.repermorder);
        if any(tmp~=T.size)
            print(['T.subsize.other and T.subsize.segment do not '...
                   'correspond with T.size.\n']);
            return;
        end
    end
    % size and repermorder do not agree
    if T.ispermuted
        if any(T.repermorder~=1:numel(T.repermorder))
            print('If T.ispermuted is true, T.repermorder should equal 1:K.\n')
            return;
        end
    end
  case 'decimate'
    
    % fields missing
    if ~isfield(T, 'size') || ~isfield(T, 'subsize') || ...
            ~isfield(T,'repermorder') || ~isfield(T,'ispermuted') || ...
            ~isfield(T,'val') || ~isfield(T,'dim') || ...
            ~isfield(T,'order') || ~isfield(T,'nsamples') || ...
            ~isfield(T,'subsample') || ~isfield(T,'shift')
        print(['T should have at least the fields val, dim, order, ',...
               'ind, size, subsize, repermorder, ispermuted, nsamples, subsample and shift\n']);
        return;
    end
    % no other/segment subsize field
    if ~isfield(T.subsize,'segment') || ~isfield(T.subsize,'other')
        print(['T.subsize should have at least the fields segment and ',...
               'other.\n']);
        return;
    end
    % tensor for H.val
    if ~ismatrix(T.val)
        print('T.val should be a vector or a matrix.\n\')
        return;
    end
    % order < 2
    if T.order < 2
        print('T.order must be larger than or equal to 2.\n');
        return;
    end
    % number of subsample elements != order-1
    if numel(T.subsample)~=T.order-1
        print('The number of T.ind is not equal to T.order-1.\n');
        return;
    end
    % number of sample numbers != 1
    if numel(T.nsamples)~=1
        print('The number of segment numbers is not equal to 1.\n');
        return;
    end
    % number of shifts != T.order-1
    if numel(T.shift)~=T.order-1
        print('The number of shifts is not equal to T.order-1.\n'); 
        return;
    end
    % negative shifts
    if any(T.shift<=0)
        print('The shifts cannot be zero or negative.\n');
        return;
    end
    % number of size is smaller than order, or than order + 1
    if size(T.val,2)==1
        if numel(T.size) ~= T.order
            print(['If T.val is a vector, the dimension of T should '...
                   'equal T.order.\n']);
            return;
        end
    else
        if numel(T.size) <= T.order
            print(['If T.val is a matrix, the dimension of T should '...
                   'exceed T.order.\n']);
            return;
        end
    end
    % number of subsize.segment elements is not equal to order
    if numel(T.subsize.segment)~=T.order
        print('The number of elements in T.subsize.segment should equal T.order.\n');
        return;
    end
    % prod(T.subsize.other) should equal size(T.val,2)
    if ~isempty(T.subsize.other) && prod(T.subsize.other)~=size(T.val,2)
        print('prod(T.subsize.other) is not equal to size(T.val,2).\n'); 
        return;
    end
    % sum(T.subsize.segment)-order+1 should equal N
    if sum(T.subsize.segment)-T.order+1 ~= size(T.val,1)
        print('sum(T.subsize.segment)-T.order+1 should equal N.\n');
        return; 
    end
    % concatenation of subsize.segment and subsize.other, and then
    % permute, is equal to size
    if T.ispermuted
        if any([T.subsize.segment T.subsize.other]~=T.size)
            print(['T.subsize.other and T.subsize.segment do not '...
                   'correspond with T.size.\n']);
            return;
        end
    else
        tmp = [T.subsize.segment T.subsize.other];
        tmp = tmp(T.repermorder);
        if any(tmp~=T.size)
            print(['T.subsize.other and T.subsize.segment do not '...
                   'correspond with T.size.\n']);
            return;
        end
    end
    % size and repermorder do not agree
    if T.ispermuted
        if any(T.repermorder~=1:numel(T.repermorder))
            print('If T.ispermuted is true, T.repermorder should equal 1:K.\n')
            return;
        end
    end
  otherwise
    print('Not a known tensor type');
    return;
end

isvalid = true;

function print(varargin)
    if printerrors
        fprintf(varargin{:});
    end
end

end
