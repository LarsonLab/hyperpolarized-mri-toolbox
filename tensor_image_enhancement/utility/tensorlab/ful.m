function T = ful(T,varargin)
%FUL Convert formatted or structured data set to an array.
%   T = ful(T) converts the formatted or structured data set T into a MATLAB
%   array. If T is incomplete, unknown entries are represented by the value NaN.
%
%   T = ful(T,ind) computes only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = ful(T,i,j,k,...) computes only the elements T(i,j,k,...) from the
%   tensor. The number of indices should match the order of the tensor, given
%   by getorder(T). The colon operator can be used as a string, e.g., for a
%   third order tensor, ful(T,i,':',k) can be used.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

type = getstructure(T);
switch type
    case 'full',
        if nargin > 1
            T = T(varargin{:});
        end
    case 'cpd',   T = cpdgen(T,varargin{:});
    case 'lmlra', T = lmlragen(T{1},T{2},varargin{:});
    case 'btd',   T = btdgen(T,varargin{:});
    case 'tt',    T = ttgen(T,varargin{:});
    case 'hankel'
        if nargin==1
            H = hankelize(T.val,'dim',1,'order',T.order,'ind',T.ind,...
                'full',true,'perm',true,'fullLimit',Inf);
            H = reshape(H,[T.subsize.hankel T.subsize.other]);
            T = permute(H,T.repermorder);
        elseif nargin==2
            if ischar(varargin{1}) && strcmp(varargin{1}, ':')
                T = reshape(ful(T),[],1);
                return;
            end
            sub = cell(1,getorder(T));
            [sub{:}] = ind2sub(getsize(T),varargin{1}(:));
        elseif nargin == getorder(T)+1 % subscripts
            for n = 1:length(varargin)
                if ischar(varargin{n}) && strcmp(varargin{n}, ':')
                    varargin{n} = 1:T.size(n);
                end
            end
            sub = cell(1,length(varargin));
            [sub{:}] = ndgrid(varargin{:});
            sub = cellfun(@(x) x(:),sub,'UniformOutput',false);
        else
            error('ful:index', ...
                'Either linear or subscripts indices should be provided');
        end
        if nargin>1
            otherdims = [1:T.dim-1 T.dim+T.order:getorder(T)];
            dims = T.dim:T.dim+T.order-1;
            dataind = sum(cat(2,sub{dims}),2)-T.order+1;
            if numel(otherdims)==0
                dataotherind = ones(size(dataind));
            elseif numel(otherdims)==1
                dataotherind = sub{otherdims};
            else
                dataotherind = sub2ind(T.subsize.other,sub{otherdims});
            end
            
            valind = sub2ind(size(T.val),dataind,dataotherind);
            
            if nargin==2
                T = reshape(T.val(valind),size(varargin{1}));
            else
                T = reshape(T.val(valind),cellfun(@numel,varargin));
            end
        end
    case 'loewner'
        if nargin==1
            L = loewnerize(T.val,'dim',1,'order',T.order,'ind',T.ind,...
                'full',true,'perm',true,'t',T.t,'fullLimit',Inf);
            L = reshape(L,[T.subsize.loewner T.subsize.other]);
            T = permute(L,T.repermorder);
        elseif nargin==2
            if ischar(varargin{1}) && strcmp(varargin{1}, ':')
                T = reshape(ful(T),[],1);
                return;
            end
            sub = cell(1,getorder(T));
            [sub{:}] = ind2sub(getsize(T),varargin{1}(:));
        elseif nargin == getorder(T)+1 % subscripts
            for n = 1:length(varargin)
                if ischar(varargin{n}) && strcmp(varargin{n}, ':')
                    varargin{n} = 1:T.size(n);
                end
            end
            sub = cell(1,length(varargin));
            [sub{:}] = ndgrid(varargin{:});
            sub = cellfun(@(x) x(:),sub,'UniformOutput',false);
        else
            error('ful:index', ...
                'Either linear or subscripts indices should be provided');
        end
        if nargin>1
            otherdims = [1:T.dim-1 T.dim+T.order:getorder(T)];
            dims = T.dim:T.dim+T.order-1;
            
            if numel(otherdims)==1
                dataotherind = sub{otherdims};
            elseif numel(otherdims)>1
                dataotherind = sub2ind(T.subsize.other,sub{otherdims});
            end
            
            totaldata = 0;
            for i = 1:numel(dims)
                idx = T.ind{i}(sub{dims(i)}).';
                
                if numel(otherdims)==0, linidx = idx;
                else linidx = sub2ind(size(T.val),idx,dataotherind); end
                
                data = T.val(linidx);
                for o = [1:i-1 i+1:numel(dims)]
                    data = data./(T.t(T.ind{i}(sub{dims(i)}))-T.t(T.ind{o}(sub{dims(o)})));
                end
                totaldata = totaldata + data;
            end
            
            if nargin==2
                T = reshape(totaldata,size(varargin{1}));
            else
                T = reshape(totaldata,cellfun(@numel,varargin));
            end
        end
        
    case {'incomplete','sparse'}
        val = T.val;
        size_tens = T.size;
        if nargin == 1
            if 8*prod(size_tens) > 8e9
                error('ful:size','T is too large to be stored as an array.');
            end
            if ~isfield(T,'ind')
                ind = sub2ind(size_tens,T.sub{:});
            else
                ind = T.ind;
            end
            if strcmp(getstructure(T),'sparse'), T = zeros(size_tens);
            else T = nan(size_tens);
            end
            T(ind) = val;
        elseif nargin == 2 % linear indices
            if ischar(varargin{1}) && strcmp(varargin{1}, ':')
                T = reshape(ful(T),[],1);
                return;
            end
            if ~isfield(T, 'ind'), T = fmt(T); end
            if strcmp(type, 'incomplete')
                val = nan(size(varargin{1}));
            else
                val = zeros(size(varargin{1}));
            end
            [ia,ib] = ismember(varargin{1}(:), T.ind);
            val(ia) = T.val(ib(ia));
            T = val;
        elseif nargin == length(size_tens)+1 % subscripts
            for n = 1:length(varargin)
                if ischar(varargin{n}) && strcmp(varargin{n}, ':')
                    varargin{n} = 1:T.size(n);
                end
            end
            size_result = cellfun(@length, varargin);
            if 8*prod(size_result) > 8e9
                error('ful:size','T is too large to be stored as an array.');
            end
            if ~isfield(T, 'ind'), T = fmt(T); end
            if strcmp(type, 'incomplete'), val = nan(size_result);
            else val = zeros(size_result); end
            sub = cell(1,length(varargin));
            [sub{:}] = ndgrid(varargin{:});
            ind = sub2ind(T.size, sub{:});
            [ia,ib] = ismember(ind(:), T.ind);
            val(ia) = T.val(ib(ia));
            T = val;
        else
            error('ful:index', ...
                'Either linear or subscripts indices should be provided');
        end
    case 'segment'
        X = desegmentize(T);
        X = segmentize(X,'dim',T.dim,'order',T.order,'segsize',T.segsize,...
            'perm',T.ispermuted,'shift',T.shift);
        T = ful(X,varargin{:});
        
    case 'decimate'
        X = dedecimate(T);
        X = decimate(X,'dim',T.dim,'order',T.order,'subsample',T.subsample,...
            'perm',T.ispermuted,'shift',T.shift);
        T = ful(X,varargin{:});
    otherwise
        error('ful:T', 'Unknown structured tensor type.');
end

