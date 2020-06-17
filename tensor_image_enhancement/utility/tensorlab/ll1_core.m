function [U,output] = ll1_core(T,U0,varargin)
%LL1_CORE computational routines for LL1 decomposition
%   LL1_CORE should not be called directly. Use LL1_MINF or LL1_NLS instead. 
%   
%   See also: ll1_minf, ll1_nls

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

    type = getstructure(T);
    unstructuredtypes = {'full', 'incomplete', 'sparse'};
    isstructured = ~any(cellfun(@(s) strcmpi(type, s), ...
                                unstructuredtypes));
    if ~isstructured, 
        T = fmt(T,true); 
        type = getstructure(T);
    end
    size_tens = getsize(T);
    sz = getsize(U0);
    size_tens = [size_tens ones(1, length(sz)-length(size_tens))];
    
    isincomplete = strcmpi(type,'incomplete');

    if ~isempty(varargin) && isnumeric(varargin{1})
        L = varargin{1};
        varargin = varargin(2:end);
    else 
        L = [];
    end
    L = L(:).';
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('OutputFormat', 'auto');
    p.addOptional('LargeScale', 'auto');
    p.addOptional('Algorithm', @nls_gndl); % or @nls_gncgs,@nls_lm
    p.addOptional('OptimizationType', 'nls'); % or @nls_gncgs,@nls_lm
    p.addOptional('CGMaxIter', 10);
    p.addOptional('Display', 0);
    p.addOptional('TolLargeScale', 0.02);
    p.addOptional('TolAbs', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    
    % Set absolute tolerance
    cache = struct;
    try 
        cache.T2 = frob(T,'squared');
    catch e
        if strcmpi(e.identifier, 'frob:notImplemented');
            error('ll1_core:notImplemented', ...
                  ['ll1_core does not support the structured tensor type %s, yet. Use ' ...
                   'ful(T) instead.'], type);
        end
    end
    if any(strcmpi(p.UsingDefaults, 'TolAbs')) 
        if isstructured
            options.TolAbs = 0.5*1e-15*cache.T2;
        else 
            options.TolAbs = 0.5*options.TolAbs*cache.T2;
        end
    end

    
    % Only third order
    if length(size_tens) ~= 3
        error('ll1_core:T', 'Only third-order tensors can be handled');
    end
    
    % convert to internal format
    if all(cellfun(@isnumeric, U0)) 
        % CPD format
        U0 = U0(:).';
        if isempty(L)
            error('ll1_core:noL', ...
                  'L should be given if U0 is in the CPD format');
        end
        if ~all(cellfun('size', U0, 2) == [sum(L), sum(L), length(L)])
            error('ll1_core:U0', ['size(U0{n},2) should be sum(L) for n=1,2 ' ...
                                'and length(L) for n=3.']);
        end
        inputformat = 'cpd';
    elseif all(cellfun(@iscell, U0)) && ... 
            all(cellfun(@(u) all(cellfun(@isnumeric,u)), U0)) 
        % BTD format: convert to cpd format
        U0 = U0(:).';
        U0 = cellfun(@(u) u(:).', U0, 'UniformOutput', false);
        U = cell(1,3);
        for k = 1:3, U{k} = cell(1, length(U0)); end
        for r = 1:length(U0)
            if ~ismatrix(U0{r}{end}) || ...
                    ~all(cellfun('size',U0{r}(1:2),2)==size(U0{r}{end})) 
                error('ll1_core:U0', ['size(U0{r}{n},2) should be size(U0{r}{end},n) ' ...
                                    'for n=1,2,3 and U0{r}{end} should be a ' ...
                                    'matrix, for all r=1:length(U0).']);
            end
            U{1}{r} = U0{r}{1}*U0{r}{end};
            U{2}{r} = U0{r}{2};
            U{3}{r} = U0{r}{3};
        end
        for k = 1:3
            if ~all(cellfun('size',U{k},1) == size_tens(k))
                error('ll1_core:U0', ['size(U0{r}{n},1) should be size(T,n) ' ...
                                    'for all n=1:3 and r=1:length(U0).']);
            end
        end
        if isempty(L)
            L = cellfun('size', U{1}, 2);
        elseif length(U{1}) ~= length(L) || ~all(cellfun('size', U{1},2)==L)
            error('ll1_core:L', ['The detected L=cellfun(@(u) size(u{end},1), U0) ' ...
                                'and the given L are not equal. (L is optional ' ...
                                'in the BTD format.)']);
        end
        U0 = cellfun(@(u) cat(2, u{:}), U, 'UniformOutput', false);
        inputformat = 'btd';
    else 
        error('ll1_core:U0', 'Unknown initialization format');
    end
    
    switch options.OutputFormat,
      case 'auto'
        options.OutputFormat = inputformat;
      case {'cpd', 'btd'}
        % ok
      otherwise
        error('ll1_core:output_format', ['The output format should be ''auto'', ' ...
                            '''cpd'' or ''btd''']);
    end 
    
    R = length(L);
    sL = sum(L);
    N = 3;
    
    if ischar(options.LargeScale)
        options.LargeScale = sum(size_tens)*sL > 100;
    end
    
    % Select solver
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
            error('ll1:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the nls optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), nlsfun)
            warning('ll1:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      case 'minf'
        if ~isfield(options, 'Algorithm')
            options.Algorithm = @minf_lbfgsdl;
        elseif ismember(func2str(options.Algorithm), nlsfun)
            error('ll1:incompatibleAlgorithm', ...
                  ['The %s method is incompatible with the minf optimization ' ...
                   'type.'], func2str(options.Algorithm));
        elseif ~ismember(func2str(options.Algorithm), minffun)
            warning('ll1:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
                                'at your own risk'], options.Algorithm);
        end
      otherwise
        error('ll1:unknownOptimizationType', ...
              ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
               'are supported now'], options.OptimizationType);
    end

    % Select optimization subroutines.
    usestate = false;
    if strcmpi(options.OptimizationType, 'nls')
        usestate = true;
        if options.LargeScale, dF.JHJx = @JHJx; else dF.JHJ = @JHJ; end
        dF.JHF = @grad;
        if ~isfield(options, 'M'), options.M = 'block-Jacobi'; end
        switch options.M
          case 'block-Jacobi', dF.M = @M_blockJacobi;
          otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
        end
    else 
        dF = @grad;
    end
    state(U0,true);
    
    [U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
    output.Name = func2str(options.Algorithm);
    
    if strcmpi(options.OutputFormat, 'btd')
        U{1} = mat2cell(U{1}, size(U{1},1), L);
        U{2} = mat2cell(U{2}, size(U{2},1), L);
        U{3} = num2cell(U{3}, 1);
        Uout = cell(1,R);
        for r = 1:R
            Uout{r} = {U{1}{r}, U{2}{r}, U{3}{r}, eye(L(r))};
        end
        U = Uout;
    end
    
    % For structured types, check the relative error
    if output.info == 4 && isstructured
        warning('ll1_core:accuracy', ...
                ['Maximal numerical accuracy for structured tensors reached. The ' ...
                 'result may be improved using ful(T) instead of T.']);
    end
    
    function state(z,firstrun)
        
        if nargin == 2 && firstrun
            % Store the fraction of known elements.
            if isincomplete
                cache.scale = length(T.val)./prod(T.size);
            end
            cache.offset = cumsum([1 cellfun(@numel,z)]);
            expvec = arrayfun(@(n) ones(1, L(n))*n, 1:R, 'UniformOutput', false); 
            cache.expvec = cat(2, expvec{:});
            cache.P = sparse(1:sum(L), cache.expvec, 1);
            cache.offset = cellfun(@numel, U0);
            cache.offset = cumsum([1 cache.offset]);
            if ~usestate, return; end
        end
        
        z{N} = z{N}(:,cache.expvec);
        cache.UHU = cellfun(@(u) conj(u'*u), z, 'UniformOutput', false);

        % Optionally cache some results for the block-Jacobi preconditioner.
        if ischar(options.M) || isa(options.M,'function_handle')
            UHU = cache.UHU;
            cache.invUHU{1} = inv(UHU{2}.*UHU{3});
            cache.invUHU{2} = inv(UHU{1}.*UHU{3});
            cache.invUHU{3} = contract(inv(UHU{1}.*UHU{2}), false);
        end
    end

    function fval = objfun(z)
        % LL1 objective function.
        z{N} = z{N}(:, cache.expvec);
        if isstructured
            fval = abs(0.5*cache.T2 - real(inprod(T, z)) + 0.5*frob(z,'squared'));
        else 
            cache.residual = cpdres(T,z,'Format',false,'Type',type);
            if isstruct(cache.residual), fval = cache.residual.val;
            else fval = cache.residual; end
            fval = 0.5*(fval(:)'*fval(:));
        end
    end

    function grad = grad(z)
        if usestate, state(z); end
        if ~isstructured, E = cache.residual; end
        offset = cache.offset;
        grad = nan(offset(end)-1, 1);
        z{N} = z{N}(:, cache.expvec);
        for n = 1:N
            if isstructured,
                tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            else 
                tmp = full(mtkrprod(E, z, n));
            end
            if n == N, tmp = tmp*cache.P; end
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end 
    end

    function jhj = JHJ(z)
        UHU = cellfun(@conj, cache.UHU, 'UniformOutput', false);
        offset = cache.offset;
        expvec = cache.expvec;
        P = cache.P;
        
        jhj = zeros(offset(end)-1);
        z{N} = z{N}(:,expvec);
        
        % Off diagonal blocks
        for n = 1:N
            idxw = true(1,N);
            idxw(n) = false;
            idxn = offset(n):offset(n+1)-1;
            Wn = prod(cat(3, UHU{idxw}),3);
            if n == N, Wn = contract(Wn, true); end
            jhj(idxn,idxn) = kron(Wn, eye(size_tens(n)));
            for m = n+1:N
                idxw = true(1,N);
                idxw([n m]) = false;
                idxn = offset(n):offset(n+1)-1;
                idxm = offset(m):offset(m+1)-1;
                Wnm = reshape(UHU{idxw}, [1 sL 1 sL]);
                JHJnm = bsxfun(@times, reshape(z{n}, [size_tens(n) 1 1 sL]), ...
                               reshape(z{m}', [1 sL size_tens(m) 1]));
                JHJnm = bsxfun(@times, JHJnm, Wnm);
                if m < N, 
                    JHJnm = reshape(JHJnm,[size_tens(n)*sL size_tens(m)*sL]);
                else 
                    JHJnm = reshape(JHJnm,[size_tens(n)*sL*size_tens(m),sL]);
                    JHJnm = JHJnm*P;
                    JHJnm = reshape(JHJnm,[size_tens(n)*sL size_tens(m)*R]);
                end
                jhj(idxn,idxm) = JHJnm;
                jhj(idxm,idxn) = JHJnm';
            end 
        end
        
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            JHJ = JHJ*cache.scale;
        end
    end
    
    function y = JHJx(z,x)
        
        offset = cache.offset;
        expvec = cache.expvec;
        UHU = cache.UHU;
        XHU = cell(1,N);
        y = nan(offset(end)-1,1);

        % Diagonal blocks
        for n = 1:3
            idxw = true(1,N);
            idxw(n) = false;
            Wn = prod(cat(3,UHU{idxw}),3);
            if n == 3, Wn = contract(Wn, true); end
            tmp = reshape(x(offset(n):offset(n+1)-1), size_tens(n), []);
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*Wn;
        end
        XHU{3} = XHU{3}(expvec,expvec);
        
        % Off diagonal blocks
        for n = 1:N-1
            idxn = offset(n):offset(n+1)-1;
            Wn = zeros(sL);
            for m = n+1:N
                idxm = offset(m):offset(m+1)-1;
                idxw = true(1,N);
                idxw([n m]) = false;
                Wn = Wn + UHU{idxw}.*XHU{m};
                Wnm = UHU{idxw}.*XHU{n};
                if m == 3, Wnm = contract(Wnm,true); end
                tmp = z{m}*Wnm;
                y(idxm) = y(idxm) + tmp(:);
            end
            tmp = z{n}*Wn;
            y(idxn) = y(idxn) + tmp(:);
        end
        
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            y = y*cache.scale;
        end
    end
    
    function x = M_blockJacobi(~,b)
    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
        x = nan(size(b));
        offset = cache.offset;
        invUHU = cache.invUHU;
        for n = 1:N
            idx = offset(n):offset(n+1)-1;
            x(idx) = reshape(b(idx), size_tens(n), [])*invUHU{n};
        end 
        
        % If incomplete, approximate the effect of missing entries.
        if isincomplete
            x = x/cache.scale;
        end
    end
    
    function W = contract(W, transposed)
        P = cache.P;
        if transposed, W = P.'*W*P;
        else W = P.'*W.'*P; end
    end

end
