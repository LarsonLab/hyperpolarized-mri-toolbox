function [U,output] = ccpd_core(dataset, varargin)
%CCPD_CORE Computational routines for coupled/symmetric CPD decomposition.
%   CCPD_CORE should not be called directly. Use CCPD_MINF or CCPD_NLS
%   instead.
%
%   See also: ccpd_minf, ccpd_nls

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/10   NV      Initial version

    if isstruct(dataset) && isfield(dataset, 'factorizations')
        % use SDF syntax
        model = sdf_check(dataset, 'internal', 'onlyCPD', ...
                          'onlySimpleFactors', 'noCPDI');
        dataset = cellfun(@(f) f.data, model.factorizations, 'UniformOutput', ...
                          false);
        cidx = cellfun(@(f) cat(2,f.factors{:}), model.factorizations, 'UniformOutput', ...
                       false);
        U0 = model.factors(unique(cat(2,cidx{:})));
        U0 = cellfun(@(u) u{1}{1}, U0, 'UniformOutput', false);
        constants = [model.isconstant{unique(cat(2, cidx{:}))}];
        for fidx = 1:length(U0)
            if ~constants(fidx)
                U0{fidx} = model.variables{U0{fidx}};
            end
        end
    elseif iscell(dataset) || isnumeric(dataset) || ...
            (isstruct(dataset) && ~isfield(dataset, 'factorizations'))
        % use old syntax
        if isnumeric(dataset), dataset = {dataset}; end
        if nargin < 2
            error('ccpd_core:U0', ['No initial guess for the factor matrices ' ...
                                'U0 provided.']);
        elseif nargin < 3
            error('ccpd_core:cidx', 'No coupling indices cidx provided');
        end
        U0 = varargin{1};
        cidx = varargin{2};
        varargin = varargin(3:end);
        if length(dataset) == 1 && ~iscell(cidx), cidx = {cidx}; end
        if ~iscell(cidx)
            error('ccpd_core:cidx', 'cidx should a cell of coupling indices');
        end
        if length(cidx) ~= length(dataset)
            error('ccpd_core:cidx', ['The cell cidx should contain one array ' ...
                                'of coupling indices for each dataset.']);
        end
        if ~iscell(U0) || any(~cellfun(@(u) isnumeric(u) && ismatrix(u), U0))
            error('ccpd_core:U0', 'U0 should be a cell of factor matrices');
        end
        cidx = cellfun(@(c) c(:).', cidx, 'UniformOutput', false);
        U0 = U0(:).';
        if any(cat(2,cidx{:}) < 1 | cat(2,cidx{:}) > length(U0))
            error('ccpd_core:cidx', ['For all n, all coupling indices i in ' ...
                                'cidx{n} should be 1 <= i <= length(U0)']);
        end

        % check consistencies with data
        size_U    = cellfun('size', U0, 1);
        size_data = cellfun(@getsize, dataset, 'UniformOutput', false);
        for f = 1:length(dataset)
            sz = size_U(cidx{f});
            szd = [size_data{f} ones(1, length(sz)-length(size_data{f}))];
            invalid = length(sz) == length(szd) && any(sz ~= szd);
            invalid = invalid || (length(sz) > length(szd) && ...
                    (any(sz ~= szd) || any(sz(length(szd)+1:end) ~= 1)));
            invalid = invalid || (length(sz) < length(szd));
            if invalid
                error('ccpd_core:dimensionMismatch', ['For all n, getsize(dataset{n}) ' ...
                                    'should match cellfun(''size'', U0(cidx{n}),1).']);
            end
        end
        if any(cellfun('size', U0, 2) ~= size(U0{1},2))
            error('ccpd_core:U0', 'All factor matrices U0 should have R columns.');
        end
    end
    
    % compute default weights
    weights = zeros(1, length(dataset));
    for f = 1:length(dataset)
        if strcmp(getstructure(dataset{f}), 'incomplete')
            weights(f) = 1/length(dataset{f}.val);
        else 
            weights(f) = 1/prod(getsize(dataset{f}));
        end
    end
    
    % fix unused factors
    idxunused = 1:length(U0);
    idxunused(cat(2,cidx{:})) = [];
    Uunused = U0(idxunused);
    [cidxorig, ~, cidxnew] = unique(cat(2,cidx{:}));
    U = U0(cidxorig);
    cidx = mat2cell(cidxnew(:).', 1, cellfun(@length, cidx));

    % set some parameters
    N = length(U);
    R = size(U{1},2);
    size_tens = cellfun('size', U, 1);
    
    % process options
    p = inputParser;
    p.addOptional('LargeScale', sum(cellfun(@numel, U)) > 1e2);
    p.addOptional('OptimizationType', 'nls');
    p.addOptional('RelWeights', ones(1, length(dataset)));
    p.addOptional('Weights', []);
    p.addOptional('Display', 0);
    p.addOptional('Debug', false);
    p.addOptional('DetectSymmetry', true);
    p.addOptional('IsSymmetric', nan(1, length(dataset)));
    p.addOptional('PC', true);
    p.addOptional('M', nan);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    
    if ~any(strcmpi(p.UsingDefaults, 'RelWeights')) && ~any(strcmpi(p.UsingDefaults, 'Weights'))
        warning('ccpd:weights',['Both relative and absolute weights ' ...
                            'are supplied, proceeding with absolute weights.']);
    end
    if isempty(options.Weights) && any(strcmpi(p.UsingDefaults, ...
                                               'RelWeights'))
        % check model for weights
        if exist('model', 'var')
            w = nan(1, length(model.factorizations));
            rw = nan(1, length(model.factorizations));
            for f = 1:length(model.factorizations)
                if isfield(model.factorizations{f}, 'weight')
                    w(f) = model.factorizations{f}.weight;
                elseif isfield(model.factorizations{f}, 'relweight')
                    rw(f) = model.factorizations{f}.relweight;
                end
            end
            if all(~isnan(w))
                options.Weights = w(:).';
            elseif all(~isnan(rw))
                options.RelWeights = rw(:).';
            end 
        end
    end
    if isempty(options.Weights)
        if length(options.RelWeights) ~= length(dataset)
            error('ccpd:weights',['The number of weights must equal the ' ...
                                'number of factorizations.']);
        end
        weights = 2*options.RelWeights(:).'/sum(options.RelWeights(:)).*weights(:).';
    else 
        weights = options.Weights(:).';
        if length(options.Weights) ~= length(dataset)
            error('ccpd:weights',['The number of weights must equal the ' ...
                                'number of factorizations.']);
        end
    end 
        
    cache = struct();
    state(U,true);
    
    if ischar(options.Debug)
        switch options.Debug
          case 'grad' 
            objfun(U0);
            U = grad(U0);
            output = [];
            return;
          case 'jhj'
            U = JHJ(U0);
            output = [];
            return;
          case 'jhjx'
            x = randn(cache.offset(end)-1,1);
            if any(cellfun(@(u) any(~isreal(u(:))), U0))
                x = x + randn(cache.offset(end)-1,1)*1i;
            end
            U = JHJx(U0,x);
            output = x;
            return;
        end
    end

    options.OptimizationType = lower(options.OptimizationType);
    nlsfun = cellfun(@func2str, {@nls_gndl,@nls_gncgs,@nls_lm}, 'UniformOutput', ...
                     false);
    minffun = cellfun(@func2str, {@minf_lbfgsdl,@minf_lbfgs,@minf_ncg}, ...
                      'UniformOutput', false);
    switch options.OptimizationType
      case 'nls'
        if ~isfield(options, 'Algorithm')
            options.Algorithm = @nls_gndl;
        elseif ismember(func2str(options.Algorithm), minffun);
            error('ccpd:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the nls optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), nlsfun)
            warning('ccpd:unknownAlgorithm', ['The %s is not known by SDF. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      case 'minf'
        if ~isfield(options, 'Algorithm')
            options.Algorithm = @minf_lbfgsdl;
        elseif ismember(func2str(options.Algorithm), nlsfun)
            error('ccpd:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the minf optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), minffun)
            warning('ccpd:unknownAlgorithm', ['The %s is not known by SDF. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      otherwise
        error('ccpd:unknownOptimizationType', ...
              ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
               'are supported now'], options.OptimizationType);
    end

    %% Run optimization algorithm.
    if strcmpi(options.OptimizationType, 'nls')
        dF.JHF = @grad;
        
        if options.LargeScale, dF.JHJx = @JHJx;
        else dF.JHJ = @JHJ; end
        if isfield(options, 'M') && ischar(options.M) 
            switch options.M
              case 'block-Jacobi', options.M = @M_blockJacobi; 
            end
        end
        if options.PC && isfield(options,'M') && isa(options.M,'function_handle')
            dF.M = options.M;
        end 
    else 
        dF = @grad;
        if any(isnan(options.M)) || ischar(options.M)
            options = rmfield(options, 'M'); 
        end
    end
    [U,output] = options.Algorithm(@objfun,dF,U,options);
    output.Name = func2str(options.Algorithm);
    output.Symmetry = cache.issymmetric;
    output.Weights = weights;
    
    % add unused variables again
    U([cidxorig, idxunused]) = [U, Uunused];

    function state(z, firstrun)
        if nargin > 1 && firstrun
            nonstructures = {'full', 'incomplete', 'sparse'};
            cache.isstructured = ~cellfun(@(d) any(strcmp(getstructure(d), ...
                                                          nonstructures)), dataset);
            cache.issymmetric = false(1, length(dataset));
            for k = 1:length(dataset)
                try 
                    cache.T2{k} = frob(dataset{k},'squared');
                catch e
                    if strcmpi(e.identifier, 'frob:notImplemented');
                        error('ccpd_core:notImplemented', ...
                              ['ccpd_core does not support the structured tensor type %s, yet. Use ' ...
                               'ful(T) instead in factorization %d.'], type, k);
                    end
                end
                cache.multiplier{k} = histc(cidx{k},1:length(z)); 
                cache.multiplier{k} = cache.multiplier{k}(:).';
                cache.uniqueidx{k} = find(cache.multiplier{k});
                
                % detect symmetry
                skip = false;
                if exist('model', 'var') % given in model
                    if isfield(model.factorizations{k}, 'issymmetric')
                        cache.issymmetric(k) = ...
                            model.factorizations{k}.issymmetric;
                        skip = true;
                    end
                end
                if ~isnan(options.IsSymmetric(k)) % given as option
                    cache.issymmetric(k) = options.IsSymmetric(k);
                    skip = true;
                end
                
                if strcmp(getstructure(dataset{k}), 'full') && ~skip && ...
                        options.DetectSymmetry 
                    cache.issymmetric(k) = true;
                    for n = cache.uniqueidx{k}
                        if cache.multiplier{k}(n) > 1
                            idx = 1:ndims(dataset{k});
                            idxn = find(cidx{k} == n);
                            p = perms(idxn(end:-1:1));
                            for l = 2:size(p,1)
                                idxt = idx;
                                idxt(idxn) = p(l,:);
                                tmp = dataset{k};
                                tmp = tmp - permute(tmp, idxt);
                                if frob(tmp)/sqrt(cache.T2{k}) > 1e-15 
                                    cache.issymmetric(k) = false; 
                                end
                            end
                        end
                    end
                end
            end
            cache.offset = cumsum([1 cellfun(@numel, z)]);
        end
        UHU = cellfun(@(u) u'*u, z, 'UniformOutput', false);
        UHU = cellfun(@(u) u(:), UHU, 'UniformOutput', false);
        UHU = cat(2,UHU{:});
        cache.UHU = UHU;
        
        if isfield(options, 'M') && (ischar(options.M) || isa(options.M,'function_handle'))
            cache.invW = repmat({zeros(R)},1,N);
            for k = 1:length(dataset)
                for n = 1:N
                    mu = cache.multiplier{k};
                    mu(n) = mu(n)-1;
                    tmp = bsxfun(@power, UHU, mu);
                    Wn = prod(tmp,2);
                    Wn = conj(reshape(Wn, [R,R]));
                    cache.invW{n} = cache.invW{n}+weights(k)*(mu(n)+1)*Wn;
                end
            end
            for n = 1:N
                cache.invW{n} = inv(cache.invW{n});
            end
        end
    end

    function fval = objfun(z)
        fval = 0;
        for k = 1:length(dataset)
            if cache.isstructured(k)
                tmp = 0.5*cache.T2{k} - real(inprod(dataset{k},z(cidx{k}))) + ...
                      0.5*frob(z(cidx{k}), 'squared');
                fval = fval + weights(k)*tmp;
            else 
                E = cpdres(dataset{k},z(cidx{k}));
                if isstruct(E)
                    fval = fval + 0.5*weights(k)*(E.val(:)'*E.val(:));
                else 
                    fval = fval + 0.5*weights(k)*(E(:)'*E(:));
                end
                cache.residual{k} = E;
            end
        end
    end
    
    function grad = grad(z)
        state(z);
        offset = cache.offset;
        grad = zeros(offset(end)-1,1);
        for k = 1:length(dataset)
            multiplier = cache.multiplier{k};
            [idx,i] = unique(cidx{k});
            if cache.isstructured(k), 
                T = dataset{k};
            else 
                E = cache.residual{k}; 
            end
            zk = z(cidx{k});
            if cache.issymmetric(k)
                for n = 1:length(idx)
                    idxn = offset(idx(n)):offset(idx(n)+1)-1;  
                    if cache.isstructured(k)
                        tmp = mtkrprod(zk,zk,i(n)) - mtkrprod(T,zk,i(n));
                    else 
                        tmp = full(mtkrprod(E,zk,i(n)));
                    end
                    tmp = weights(k)*multiplier(idx(n))*tmp;
                    grad(idxn) = grad(idxn) + tmp(:);
                end
            else 
                idx = cidx{k};
                for n = 1:length(cidx{k})
                    idxn = offset(idx(n)):offset(idx(n)+1)-1;  
                    if cache.isstructured(k)
                        tmp = mtkrprod(zk,zk,n) - mtkrprod(T,zk,n);
                    else 
                        tmp = full(mtkrprod(E,zk,n));
                    end
                    tmp = weights(k)*tmp;
                    grad(idxn) = grad(idxn) + tmp(:);
                end
            end
        end
    end
    
    function jhj = JHJ(z)
        UHU = cache.UHU;
        jhj = zeros(cache.offset(end)-1);
        for n = 1:length(z)
            idxn = cache.offset(n):cache.offset(n+1)-1;
            for k = 1:length(dataset)
                multiplier = cache.multiplier{k}(:)';
                freq = cache.multiplier{k}(:);
                
                if multiplier(n) == 0; continue; end 
                multiplier(n) = multiplier(n) - 1;
                tmp = bsxfun(@power, UHU, multiplier);
                Wn = freq(n)*prod(tmp,2);
                Wn = reshape(Wn, [R,R]);
                jhj(idxn,idxn) = jhj(idxn,idxn) + weights(k)*kron(Wn, eye(size_tens(n)));
                
                if multiplier(n) >= 1, 
                    % update in the case of symmetries
                    multiplier(n) = multiplier(n) - 1;
                    tmp = bsxfun(@power, UHU, multiplier);
                    Wn = prod(tmp,2);
                    Wn = reshape(Wn, [R,R]);
                    JHJnn = bsxfun(@times, reshape(z{n},[size_tens(n) 1 1 R]),...
                                   reshape(z{n}',[1 R,size_tens(n),1]));
                    JHJnn = bsxfun(@times, JHJnn, reshape(Wn,[1 R 1 R]));
                    JHJnn = weights(k)*reshape(JHJnn,[size_tens(n)*R,size_tens(n)*R]);
                    jhj(idxn,idxn) = jhj(idxn,idxn) + freq(n)*(freq(n)-1)*JHJnn;
                    multiplier(n) = multiplier(n) + 1;
                end 
                
                M = n + find(freq(n+1:end)' > 0);
                for m = M(:).'
                    idxm = cache.offset(m):cache.offset(m+1)-1;
                    multiplier(m) = multiplier(m) - 1;
                    tmp = bsxfun(@power, UHU, multiplier);
                    Wnm = freq(n)*freq(m)*prod(tmp,2);
                    JHJnm = bsxfun(@times, reshape(z{n},[size_tens(n) 1 1 R]),...
                                   reshape(z{m}',[1 R,size_tens(m),1]));
                    JHJnm = bsxfun(@times, JHJnm, reshape(Wnm,[1 R 1 R]));
                    JHJnm = weights(k)*reshape(JHJnm, [size_tens(n)*R, size_tens(m)*R]);
                    jhj(idxn,idxm) = jhj(idxn,idxm) + JHJnm;
                    jhj(idxm,idxn) = jhj(idxm,idxn) + JHJnm';
                    % undo multiplier decrease
                    multiplier(m) = multiplier(m) + 1;
                end
            end
        end
    end
    
    
    function y = JHJx(z,x)
        offset = cache.offset;
        UHU = conj(cache.UHU);
        XHU = zeros(R,R,N); 
        y = zeros(offset(end)-1,1);
        for k = 1:length(dataset)
            freq = cache.multiplier{k}(:);
            for n = 1:length(z)
                if freq(n) == 0, continue; end
                multiplier = cache.multiplier{k}(:)';
                multiplier(n) = multiplier(n) - 1;
                tmp = bsxfun(@power, UHU, multiplier);
                Wn = prod(tmp,2);
                tmp = reshape(x(offset(n):offset(n+1)-1),size_tens(n),R);
                XHU(:,:,n) = freq(n)*conj(tmp'*z{n});
                idxn = offset(n):offset(n+1)-1;
                JHJnx = weights(k)*freq(n)*(tmp*reshape(Wn,[R R]));
                y(idxn) = y(idxn) + JHJnx(:);
            
                if freq(n) > 1 % add off diagonal blocks
                    multiplier(n) = multiplier(n) - 1;
                    tmp = bsxfun(@power, UHU, multiplier);
                    Wnm = reshape(prod(tmp,2), R, R);
                    Wn = (freq(n)-1)*(Wnm.*XHU(:,:,n));
                    JHJnx = weights(k)*z{n}*Wn;
                    y(idxn) = y(idxn) + JHJnx(:);
                end
            end
            
            for n = 1:length(z)-1
                if freq(n) == 0, continue; end
                multiplier = cache.multiplier{k}(:)';
                multiplier(n) = multiplier(n) - 1;
                idxn = offset(n):offset(n+1)-1;
                Wn = zeros(R);
                M = n + find(freq(n+1:end)' > 0);
                for m = M(:).';
                    multiplier(m) = multiplier(m) - 1;
                    idxm = offset(m):offset(m+1)-1;
                    if length(cidx{k}) == 2
                        Wn = Wn+XHU(:,:,m);
                        JHJmnx = z{m}*XHU(:,:,n);
                    else
                        tmp = bsxfun(@power, UHU, multiplier);
                        Wnm = prod(tmp,2);
                        Wnm = reshape(Wnm,[R R]);
                        Wn = Wn+Wnm.*XHU(:,:,m);
                        JHJmnx = z{m}*(Wnm.*XHU(:,:,n));
                    end
                    y(idxm) = y(idxm)+freq(m)*weights(k)*JHJmnx(:);
                    multiplier(m) = multiplier(m) + 1;
                end
                JHJnx = z{n}*Wn;
                y(idxn) = y(idxn)+weights(k)*freq(n)*JHJnx(:);
            end
        end
    end

    % Solve Mx = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factor matrices.
    function x = M_blockJacobi(~,b)
        offset = cache.offset;
        x = zeros(size(b));
        for n = 1:N
            idx = offset(n):offset(n+1)-1;
            x(idx) = reshape(b(idx),size_tens(n),R)*cache.invW{n};
        end
    end

end 