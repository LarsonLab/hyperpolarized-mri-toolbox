function [U, output] = cpd_rbs(T, U, varargin)
%CPD_RBS CPD by randomized block sampling
%   [U,output] = CPD_RBS(T,U) computes factor matrices U{1}, ..., U{N} belonging
%   to a canonical polyadic decomposition of the N-th order tensor T by
%   minimizing by minimizing 0.5*frob(T-cpdgen(U))^2. The algorithm is
%   initialized with the factor matrices U0{n}. The tensor T can be full, sparse
%   or structured. Incomplete tensors are not supported. The structure output
%   returns additional information:
%
%      output.derivative  - The approximate derivative of relative step size
%      output.relchange   - Relative change of variables w.r.t. CRB
%      output.delta       - Relative step radius restriction
%      output.name        - The name of the selected algorithm.
%      output.<...>       - The output of the selected algorithm.
%    
%   CPD_RBS(T,U,'key',value,...) or CPD_RBS(T,U,options) in which options is
%   a struct('key', value, ...) can be used to pass the following options:
%   
%     Sampling options:
%   
%     - BlockSize           - Size of the sampled blocks. Default:
%                             min(2*R,R+10)
%     - Shuffle             - If false, do not permute the indices when no
%                             more full blocks are left. Default: true.
%     - FullRandom          - If true, permute the indices every time a block
%                             is sampled, instead of when no more blocks are
%                             left. Default: false.
%
%     Step restriction options:
%   
%     - Restriction         - If Restriction is a number, the relative step
%                             restriction is constant for all iterations. If
%                             Restriction is a function handle, it should
%                             accept three values as input: 
%                             - delta: the previous restriction value
%                             - k: the current iteration
%                             - output: the current output struct 
%                             The function should return a single relative
%                             step restriction value. For example:
%
%                               Restriction = @(delta, k, output) 0.95^k. 
% 
%                             Default: 0.8.
%     - AutoStep            - If true, a heuristic tries to detect whether
%                             algorithm has converged to a region around the
%                             optimum, and will then start restricting the
%                             relative step size to improve accuracy, using a
%                             restriction of the form 0.95^(k/Q) in which Q
%                             is the maximum number of blocks per dimension.
%                             Default: true if the default is used for
%                             Restriction. 
%
%     Stopping criterion options:
%
%     - TolCRB              - Tolerance for the Cramér-Rao based stopping
%                             criterion. Default: 0.5.
%     - TolX                - Tolerance for the relative step size. Default:
%                             1e-8.
%     - TolFun              - Tolerance for the relative objective function
%                             value. Default: 1e-16.
%
%     Algorithm options:
%   
%     - Algorithm           - Optimization algorithm to use, e.g.,
%                             minf_lbfgsdl, nls_gndl, nls_lm. Default: nls_gndl.
%     - AlgorithmOptions    - Options struct to be passed to algorithm.
%                             Default: struct('MaxIter', 1).    
%     - MaxIter             - Maximum number of iterations for the RBS
%                             algorithm (this is different from the one in
%                             AlgorithmOptions). Default: 1000
%     - LargeScale          - If true, the large scale implementation is
%                             used. Default: prod(blocksize)*R > 100
%     - M                   - Preconditioner to use in the case of the large
%                             scale implementation. Default: 'block-Jacobi'
%                             if an NLS type algorithm is used, otherwise the
%                             default of the Algorithm is used. 
%
%     Miscellaneous: 
% 
%     - Measure             - If measure is a function handle, every
%                             ItersPerMeasure iterations Measure(U) is called,
%                             and the result is printed if Display > 0. The
%                             field output.measure contains the results of
%                             measure for each iteration it is computed, and
%                             nan for other iterations. Measure(U) should
%                             return a single, numerical value. Default: nan.
%     - ItersPerMeasure     - Every ItersPerMeasure iterations, Measure is
%                             called, if given. Default: 4.
%     - SNR                 - If non-full tensor is given, Gaussian noise with
%                             a given SNR is added to the sample using
%                             noisy(Tsample, SNR, [], [], 1e5), if SNR is
%                             finite. Default: inf.
%                             If SNR is finite, no more entries than the
%                             number of entries in the tensor can be sampled.
%                             Otherwise, the chance of sampling the same
%                             entry with diffent noise values becomes large,
%                             and the result can too accurate because of
%                             noise averaging.
%     - CRBMethod           - Method to use in the CPD_CRB algorithm. See
%                             CPD_CRB for more information.
%
%     See Also: cpd_crb, cpd, cpdgen.
    
%   Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%              Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] N. Vervliet, L. De Lathauwer, "A randomized block sampling approach
%       to canonical polyadic decomposition of large-scale tensors," J. Sel.
%       Topics Signal Process., IEEE, Vol. 10, No. 2, pp. 284-295, 2016.
%
%   Version History:
%   - 2014/06/12   NV      Initial version
    
    type = getstructure(T);
    if strcmp(type, 'incomplete')
        error('cpd_rbs:notImplemented', ...
              'Incomplete tensors are not supported, yet');
    end

    R = size(U{1}, 2);
    size_tens = cellfun(@(u) size(u, 1), U);
    N = length(size_tens);
    
    if any(cellfun('size',U,2) ~= R)
        error('cpd_rbs:U','size(U{n},2) should be the same for all n.');
    end
    if length(U) < getorder(T)
        error('cpd_rbs:U','length(U) should be equal to getorder(T).');
    end
    if any(size_tens(1:getorder(T))~=getsize(T)) || ...
            any(size_tens(getorder(T)+1:end)~=1)
        error('cpd_rbs:U','size(U{n},1) should be equal to size(T,n).');
    end

    %% Option parsing
    p = inputParser;
    p.KeepUnmatched = true;
    p.addOptional('CGMaxIter', 15);
    p.addOptional('MaxIter', 1000);
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-6);
    p.addOptional('LargeScale', 'auto');
    p.addOptional('TolCRB', 0.5);
    p.addOptional('Display', 0);
    p.addOptional('BlockSize', ones(1, N)*min(2*R,R+10));
    p.addOptional('Measure', nan);
    p.addOptional('Shuffle', true);
    p.addOptional('ItersPerMeasure', 4);
    p.addOptional('AlgorithmOptions', struct('MaxIter', 1));
    p.addOptional('Restriction', 0.8);
    p.addOptional('AutoStep', 'auto');
    p.addOptional('SNR', inf);
    p.addOptional('FullRandom', false);    
    p.addOptional('M', nan);
    p.addOptional('CRBMethod', 'diag');
    p.parse(varargin{:});
    options = p.Results;
    
    if ~isa(options.Measure, 'function_handle')
        options.ItersPerMeasure = 0;
    end
    if ~(isscalar(options.Restriction) && isnumeric(options.Restriction)) && ...
            ~isa(options.Restriction, 'function_handle')
        error('cpd_rbs:restriction', ['Restriction should be a scalar or a ' ...
                            'function handle']);
    end
    
    %% Block size options
    blocksize = options.BlockSize(:).';
    if numel(blocksize) == 1, blocksize = blocksize*ones(1, N); end    
    if length(blocksize)~=length(U)
        error('cpd_rbs:blocksize','blocksize should have length 1 or length(U)');
    end
    if any(size_tens<blocksize)
        error('cpd_rbs:blocksize','blocksize should be <= getsize(T)');
    end
    blocks    = max(ceil(size_tens./blocksize));
    Q         = prod(blocksize);
    display   = options.Display;
        
    %% step size restriction options
    if isnumeric(options.Restriction)
        deltafun = @(delta, k, output) options.Restriction;
    else
        deltafun = options.Restriction;
    end
    if ischar(options.AutoStep) && strcmpi(options.AutoStep, 'auto')
        options.AutoStep = any(strcmpi('Restriction', p.UsingDefaults));
    end
    % Horizon for the determination of the derivative
    H = -2*blocks:0;
    
    %% Set algorithm options
    algoptions = p.Results.AlgorithmOptions;
    if ischar(options.LargeScale) && strcmp(options.LargeScale, 'auto')
        options.LargeScale = Q*R > 1e2;
    end
    
    useexternalroutines = false;
    usingALS = false;
    if strncmpi(func2str(options.Algorithm), 'cpd_', 4)
        useexternalroutines = true;
        if strncmpi(func2str(options.Algorithm), 'cpd_als', 7)
            if any(strcmpi('Restriction', p.UsingDefaults))
                usingALS = true;
                deltafun = @(delta, k, output) 1;
            end
        end
    elseif strncmp(func2str(options.Algorithm), 'minf', 4)
        % first order algorithm
        dF = @grad;
        usestate = false;
    else 
        % nls type algorithm
        dF.JHF = @grad;
        if options.LargeScale
            dF.JHJx = @JHJx;
        else 
            dF.JHJ = @JHJ;
        end 
        if isnan(options.M), options.M = 'block-Jacobi'; end
        switch options.M
          case 'block-Jacobi', dF.M = @M_blockJacobi;
          otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
        end
        usestate = true;
        algoptions.CGMaxIter = options.CGMaxIter;
    end
    if ~isnan(options.M)
        algoptions.M = options.M;
    end
    
    %% Set outputs
    output = struct;
    output.fval       = [];
    output.relfval    = [];
    output.relstep    = [];
    output.measure    = [];
    output.delta          = [];
    output.rho          = [];
    output.derivative = [];
    output.relchange  = [];
    output.info       = 0;
    output.iterations = 0;
    output.name = func2str(options.Algorithm);
    output.blocksize  = blocksize;
    
    %% Algorithm settings
    k         = 1;          % iteration counter
    ctr       = size_tens;  % reshuffle counter for each mode
    shuffled  = cell(1, N); % shuffled indices
    delta         = 0.5;        % init for restriction radius
    Uprev     = U;          % Initial value
    cache     = struct;     % cache for nls methods
    bold      = '%s';       % bold setting for printing
    relchange = nan;        % relative step change
    
    %% Main algorithm loop
    while ~output.info
        % Shuffle indicies if needed 
        ctr = ctr + blocksize;
        for n = 1:N
            if options.FullRandom || ctr(n) + blocksize(n) > size_tens(n)
                ctr(n) = 0; 
                if options.Shuffle
                    shuffled{n} = randperm(size_tens(n));
                else 
                    shuffled{n} = 1:size_tens(n);
                end
            end
        end
        
        % Select sampling indices 
        idx = arrayfun(@(n) ctr(n)+(1:blocksize(n)),1:N,'UniformOutput',0);
        idx = cellfun(@(s,i) sort(s(i)), shuffled, idx, 'UniformOutput', 0);
        
        % Select subtensor and affected variables
        if isnumeric(T), 
            Ts = T(idx{:});
        else 
            Ts = ful(T, idx{:}); 
            if isfinite(options.SNR), Ts = noisy(Ts, options.SNR); end
        end
        Us = cellfun(@(u,i) u(i,:), U, idx, 'UniformOutput', false);
        
        % determine step length
        pnorm = sqrt(sum(cellfun(@frob, Us).^2));
        
        if k > -min(H) + blocks
            % compute derivative 
            derivative = mean(log10(output.relfval(k+H)));
            derivative = derivative - mean(log10(output.relfval(k+H-blocks)));
            output.derivative = [output.derivative derivative/blocks];
        else 
            output.derivative = [output.derivative nan];
        end
        if options.AutoStep
            % automatic restriction heuristic
            if k > blocks*10  % all variables updated > 10 times (approx)
                bnd = 1e-4;
                if all(output.derivative(end-3*blocks:end) > -bnd)
                    if usingALS
                        deltafun = @(~,l,~) 0.95^((l-k)/blocks);
                    else 
                        deltafun = @(~,l,~) output.relstep(end)*0.95^((l-k)/blocks);
                    end
                    options.AutoStep = false; % disable checking
                end
            end
        end
        delta = deltafun(delta, k, output);
        
        % Set state
        state(Us);
        
        % compute update
        if usingALS
            algoptions.Delta = delta;
        else 
            algoptions.Delta = delta*pnorm;
        end
        algoptions.Display = 0;
        if useexternalroutines
            [Us, out] = options.Algorithm(Ts, Us, algoptions);
        else 
            [Us, out] = options.Algorithm(@(z)objfun(Ts,z), dF, Us, ...
                                          algoptions);
        end 

        % update original factors
        for n = 1:N, U{n}(idx{n},:) = Us{n}; end

        % Update output structure
        if k > 1, 
            out.fval = out.fval(2:end); 
        end
        output.fval = [output.fval, out.fval/Q];
        output.relfval = [output.relfval, out.fval/Q/output.fval(1)];
        output.relstep = [output.relstep, out.relstep];
        output.delta = [output.delta, delta];
        if ~isfield(out, 'rho'), out.rho = nan; end
        output.rho = [output.rho, out.rho];
        
        output.iterations = output.iterations + out.iterations;
        
        if k == options.MaxIter, output.info = 3; end
        if output.relfval(end) <= options.TolFun, output.info = 1; end
        if output.relstep(end) <= options.TolX, output.info = 2; end
        if options.ItersPerMeasure > 0 && mod(k, options.ItersPerMeasure) == 0
            output.measure = [output.measure options.Measure(U)];
        else
            output.measure = [output.measure nan];
        end
        
        if k > 2*blocks && mod(k, blocks) == 0
            % compute CRB
            fval = output.fval(end-min(length(output.fval),2*blocks):end);
            crb = cpd_crb(U,fval,Uprev,'Method',options.CRBMethod);
            
            % compute stopping criterion
            [~, P, D] = cpderr(Uprev, U);
            u1 = cellfun(@(u, d) u*P*d, U, D, 'UniformOutput', false);
            u1 = cellfun(@(u) u(:), u1, 'UniformOutput', false);
            u2 = cellfun(@(u) u(:), Uprev, 'UniformOutput', false);
            relchange = mean(abs(cat(1,u1{:})-cat(1,u2{:}))./sqrt(crb));
            if relchange < options.TolCRB && k > 10*blocks
                output.info = 4;
            end
            % update outputs 
            Uprev = U;
        end 
        
        output.relchange(end+1) = relchange;

        % call print if wanted
        if display > 0, print(false); end

        % Next iteration
        k = k + 1;
    end
    
    if display > 0 && mod(output.iterations, display) > 0; print(true); end
    
    if ~isnumeric(T) && isfinite(options.SNR) && ...
            prod(blocksize)*k > prod(size_tens)
        warning('cpd_rbs:noise', ['Results may be too accurate because of the ' ...
                            'noise addition']);
    end
      
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin auxiliary functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    function print(finished)
    % Display progress.
        if output.iterations == 1
            bold = '%s';
            [~,~,~,~,v] = regexp(version('-release'),'([0-9]+)([ab])');
            if usejava('Desktop') && str2double(v{1}{1}) > 2011 || ...
                    (str2double(v{1}{1}) == 2011 && strcmpi(v{1}{2},'b'))
                bold = '<strong>%s</strong>';
            end
        end
        if output.iterations == 1 || ...
                mod(output.iterations,15*display) == 0
            fprintf('\n%7s%s','', sprintf(bold,'fval'));
            fprintf('%13s%s', '', sprintf(bold,'relfval'));
            fprintf('%10s%s', '', sprintf(bold,'relstep'));
            fprintf('%10s%s',  '', sprintf(bold,'Delta'));
            fprintf('%8s%s', '', sprintf(bold,'Change'));
            if options.ItersPerMeasure > 0
                fprintf('%7s%s', '', sprintf(bold,'Measure'));
            end 
            fprintf('\n%21s%9s = %4.e %6s = %4.e\n\n','=1/2*norm(F)^2', ...
                    'TolFun',options.TolFun,'TolX',options.TolX);
        end
        if output.iterations == 1
            fprintf('%4i: % 14.8e |\n',0,output.fval(1));
        end
        if mod(output.iterations, display) ==  0 || finished
            if options.ItersPerMeasure > 0 && isnan(output.measure(end))
                output.measure(end) = options.Measure(U);
            end
            
            if options.ItersPerMeasure > 0
                fmtstring = '%4i: % 14.8e | %14.8e | %14.8e | %8.4e | %8.4e | %14.8e\n';
                fprintf(fmtstring, ...
                        output.iterations,output.fval(end), ...
                        output.relfval(end),output.relstep(end), ...
                        delta, ...
                        relchange, ...
                        output.measure(end));
            else 
                fmtstring = '%4i: % 14.8e | %14.8e | %14.8e | %8.4e | %8.4e\n';
                fprintf(fmtstring, ...
                        output.iterations,output.fval(end), ...
                        output.relfval(end),output.relstep(end), ...
                        delta, ...
                        relchange);            
            end
        end
    end
    
function state(z,~)

    % Cache the factor matrices' Gramians.
    cache.UHU = zeros(N,R*R);
    for n = 1:N
        tmp = conj(z{n}'*z{n});
        cache.UHU(n,:) = tmp(:);
    end
    
    cache.offset = cumsum([1 cellfun(@numel, z)]);
    
    % Optionally cache the inverses of the Gramians for the preconditioner.
    % In a faster language, this should be the Cholesky factor instead.
    if ischar(options.M) || isa(options.M,'function_handle')
        cache.invW = cell(1,N);
        for n = 1:N
            tmp = cache.UHU([1:n-1 n+1:N],:);
            if N > 2, tmp = prod(tmp,1); end
            cache.invW{n} = inv(reshape(tmp,[R R]));
        end
    end
end

function fval = objfun(Tb, z)
% CPD objective function.

    E = cpdres(Tb, z);
    cache.residual = E;
    if isstruct(E)
        fval = 0.5*abs(E.val(:)'*E.val(:));
    else 
        fval = 0.5*abs(E(:)'*E(:));
    end
end

function grad = grad(z)
    if usestate, state(z); end
    offset = cache.offset;    
    grad = nan(offset(end)-1,1);

    % CPDn scaled conjugate cogradient.
    E = cache.residual;
    for n = 1:length(z)
        tmp = full(mtkrprod(E,z,n));
        grad(offset(n):offset(n+1)-1) = tmp(:);
    end
end

function JHJ = JHJ(z)
    
    % CPD Jacobian's Gramian.
    UHU = conj(cache.UHU);
    JHJ = zeros(cache.offset(end)-1);
    for n = 1:N
        idxn = cache.offset(n):cache.offset(n+1)-1;
        Wn = reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]);
        JHJ(idxn,idxn) = kron(Wn,eye(size(z{n},1)));
        for m = n+1:N
            idxm = cache.offset(m):cache.offset(m+1)-1;
            Wnm = reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),[R R]);
            JHJnm = bsxfun(@times,reshape(z{n},[size(z{n},1) 1 1 R]), ...
                    reshape(conj(z{m}),[1 size(z{m},1) R 1]));
            JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 1 R R]));
            JHJnm = permute(JHJnm,[1 3 2 4]);
            JHJnm = reshape(JHJnm,[size(z{n},1)*R size(z{m},1)*R]);
            JHJ(idxn,idxm) = JHJnm;
            JHJ(idxm,idxn) = JHJnm';
        end
    end
end

function y = JHJx(z,x)
    
    % CPD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    offset = cache.offset;
    UHU = cache.UHU;
    XHU = zeros(R,R,N);
    y = nan(offset(end)-1,1);
    for n = 1:N
        Wn = UHU([1:n-1 n+1:N],:);
        if N > 2, Wn = prod(Wn,1); end
        tmp = reshape(x(offset(n):offset(n+1)-1),[],R);
        XHU(:,:,n) = conj(tmp'*z{n});
        y(offset(n):offset(n+1)-1) = tmp*reshape(Wn,[R R]);
    end
    for n = 1:N-1
        idxn = offset(n):offset(n+1)-1;
        Wn = zeros(R);
        for m = n+1:N
            idxm = offset(m):offset(m+1)-1;
            if N == 2
                Wn = Wn+XHU(:,:,m);
                JHJmnx = z{m}*XHU(:,:,n);
            else
                Wnm = UHU([1:n-1 n+1:m-1 m+1:N],:);
                if N > 3, Wnm = prod(Wnm,1); end
                Wnm = reshape(Wnm,[R R]);
                Wn = Wn+Wnm.*XHU(:,:,m);
                JHJmnx = z{m}*(Wnm.*XHU(:,:,n));
            end
            y(idxm) = y(idxm)+JHJmnx(:);
        end
        JHJnx = z{n}*Wn;
        y(idxn) = y(idxn)+JHJnx(:);
    end
end

function x = M_blockJacobi(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = nan(size(b));
    for n = 1:length(cache.offset)-1
        idxn = cache.offset(n):cache.offset(n+1)-1;
        tmp = reshape(b(idxn),[],size(cache.invW{1},1))*cache.invW{n};
        x(idxn) = tmp(:);
    end
end
end
