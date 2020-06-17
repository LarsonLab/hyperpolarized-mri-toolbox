function [U,output] = ll1(T,L,varargin)
%LL1 Decomposition in LL1 terms
%   [U,output] = LL1(T,L) computes R = length(L) terms U{r} corresponding to a
%   LL1 decomposition in BTD format of the third-order tensor T by minimizing
%   0.5*frob(T-ll1gen(U))^2. Each term U{r} is a cell array of N factor matrices
%   U{r}{n} with L(r) columns for n=1,2, and vector for n=3, followed by an
%   identity matrix U{r}{N+1}.
%
%   [U,output] = LL1(T,U0) allows the user to provide initial factors U0 in the
%   BTD format, to be used by the main algorithm. By providing an
%   initialization, the compression step will be skipped. Factor matrices equal
%   to the empty matrix [] will be initialized using the chosen initialization
%   method.
%
%   [U,output] = LL1(T,U0,L) allows the user to provide initial factors U0 in
%   the CPD format, to be used by the main algorithm. In the CPD format, the
%   first factors U0{1} and U0{2} have sum(L) columns, the third U0{3} has
%   R=length(L) columns. By providing an initialization, the compression step
%   will be skipped. Factor matrices equal to the empty matrix [] will be
%   initialized using the chosen initialization method. The output U is also
%   in the CPD format.
%
%   The structure output contains the output of the selected algorithms:
%
%      output.Preprocessing  - Preprocessing performed before compression.
%      output.Compression    - The output of the compression step.
%      output.Initialization - The output of the initialization step.
%      output.Algorithm      - The output of the main algorithm.
%      output.Refinement     - The output of the the refinement step after
%                              decompression (if applicable).
%
%   LL1(T,R,options), LL1(T,U0,options) and LL1(T,U0,L,options) allow the user
%   to choose which algorithms will be used in the different steps and also set
%   parameters for those algorithms:
%
%      Use options.Compression = [{'auto'}|@mlsvd_rsi,@lmlra_aca,@mlsvd] to
%      select whether or not the tensor is first compressed using the given
%      algorithm to a tensor of dimensions min(size(T),R*ones(1,N)). If one or
%      more initial factor matrices are provided compression is skipped. On
%      'auto', compression is only executed if it is expected to lead to a
%      significant reduction in computational complexity. In the case the tensor
%      is full or sparse, or if the tensor is not incomplete and
%      prod(getsize(T)) <= options.ExpandLimit (see further), mlsvd_rsi is the
%      default compression algorithm for 'auto'; in all other cases the default
%      is lmlra_aca. The structure options.CompressionOptions will be passed to
%      the chosen algorithm.
%
%      Use options.Initialization = [{'auto'}|@ll1_rnd|@ll1_gevd] to
%      choose the initialization method. The structure
%      options.InitializationOptions will be passed to the chosen
%      initialization method. On 'auto', @ll1_gevd is used when possible,
%      otherwise @ll1_rnd.
%
%      Use options.Algorithm = [@ll1_minf|{@ll1_nls}] to choose the main
%      algorithm. The structure options.AlgorithmOptions will be passed to the
%      chosen algorithm.
%
%      Use options.Refinement = [@ll1_minf|{@ll1_nls}|false] to choose the
%      algorithm used to refine the solution after decompression, if applicable.
%      By default, options.Refinement is equal to options.Algorithm. The
%      structure options.RefinementOptions will be passed to the chosen
%      refinement method. Parameters not set in options.RefinementOptions will
%      be copied from options.AlgorithmOptions where possible.
%
%   Further options are:
%
%      options.Display = false - Set to true to enable printing output
%                                information to the command line. If
%                                options.Display > 1, it is passed on to
%                                AlgorithmOptions and RefinementOptions,
%                                unless these structs define Display.
%      options.OutputFormat =  - Select the output format for U (either 'btd'
%      ['btd','cpd']             'cpd'. If not given, the same format as U0
%                                is chosen if U0 is given, otherwise the
%                                default is 'btd'.
%      options.ExpandLimit =   - Create full tensor if T is not incomplete
%      1e6                       and prod(getsize(T)) <= options.ExpandLimit
%                                for compression and/or initialization steps
%      options.Complex =       - If true, compute a complex solution even if
%      ['auto', true]            the tensor is real. If auto, a complex
%                                solution is only computed if T is complex,
%                                or if options.Initialization results in a
%                                complex initialization, e.g., if
%                                options.InitializationOptions.Imag=@randn.
%                                This option only works for
%                                options.Initialization = ll1_rnd or
%                                ll1_gevd.
%      options.ExploitStructure -Detect structure in the tensor and
%       = [false,true,{'auto'}]  exploit this for speed purposes. If true,
%                                this is applied for the decompositions of
%                                both the compressed and uncompressed
%                                tensor. If auto, this is only applied for
%                                the decomposition of the
%                                uncompressed tensor.
%
%   The following options are passed on to AlgorithmOptions and
%   RefinementOptions unless these structs define these options:
%      options.TolX            - Tolerance for step length
%      options.TolFun          - Tolerance for function value
%      options.MaxIter         - Maximum number of iterations
%      options.CGMaxIter       - Maximum number of CG iterations for inexact
%                                NLS type algorithms.
%
%   See also ll1gen, ll1_nls, ll1_minf, ll1_gevd, ll1_rnd, mlsvd_rsi,
%   lmlra_aca, mlsvd

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the tensor T.
type = getstructure(T);
isstructured = ~any(strcmp(type, {'full', 'sparse', 'incomplete'}));

if ~isstructured
    T = fmt(T);
    type = getstructure(T);
end
size_tens = getsize(T);

if ~isempty(varargin) && isnumeric(varargin{1})
    if ~iscell(L)
        % error
    else
        U = L;
        L = varargin{1};
        varargin = varargin(2:end);
    end
elseif iscell(L)
    U = L;
    if ~all(cellfun(@iscell, U))
        error('ll1:L', ['If U0 is given in the CPD format, L should be provided ' ...
            'as third argument']);
    end
    L = cellfun(@(u) size(u{1},2), U);
end
L = L(:).';

R = length(L);
N = length(size_tens);

if isstructured
    try 
        frob(T);
    catch e
        if strcmpi(e.identifier, 'frob:notImplemented');
            error('ll1:notImplemented', ...
                  ['ll1 does not support the structured tensor type %s, yet. Use ' ...
                   'ful(T) instead.'], type);
        end
    end
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
p = inputParser;
p.addOptional('Compression', 'auto');
p.addOptional('CompressionOptions', struct());
p.addOptional('Initialization', 'auto');
p.addOptional('InitializationOptions', struct());
p.addOptional('Algorithm', @ll1_nls); %@ll1_nls,@ll1_minf
p.addOptional('AlgorithmOptions', struct());
p.addOptional('Refinement', 'auto');
p.addOptional('RefinementOptions', struct());
p.addOptional('Display', false);
p.addOptional('OutputFormat', 'auto');
p.addOptional('Complex', 'auto');
p.addOptional('ExpandLimit', 1e6);
p.addOptional('TolX', nan);
p.addOptional('TolFun', nan);
p.addOptional('MaxIter', nan);
p.addOptional('CGMaxIter', nan);
p.addOptional('ExploitStructure','auto');
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

inputFormat = 'na';
if exist('U', 'var')
    if all(cellfun(@isnumeric, U))
        inputFormat = 'cpd';
        % do input checks
        U = U(:).';
        if length(U) ~= 3
            error('ll1:U0', 'U0 should have exactly three factors');
        end
        if ~all(cellfun(@ismatrix, U))
            error('ll1:U0', 'All entries in U0 should be matrices');
        end
        if size(U{1},2)~=sum(L) || size(U{2},2)~=sum(L) || ...
                size(U{3},2)~=length(L)
            error('ll1:U0', ['U0{1} and U0{2} should have sum(L) columns, U0{3} ' ...
                'should have length(L) columns.']);
        end
        if any(cellfun('size', U, 1)~=size_tens)
            error('ll1:U0', 'size(U0{n},1) should be size(T,n) for n=1:3');
        end
        
        % convert to btd format
        tmp = U;
        tmp{1} = mat2cell(tmp{1}, size(tmp{1},1), L);
        tmp{2} = mat2cell(tmp{2}, size(tmp{2},1), L);
        tmp{3} = num2cell(tmp{3}, 1);
        U = cell(1, R);
        for r = 1:R
            U{r} = cell(1, 4);
            U{r}{1} = tmp{1}{r};
            U{r}{2} = tmp{2}{r};
            U{r}{3} = tmp{3}{r};
            U{r}{4} = eye(L(r));
        end
    else
        inputFormat = 'btd';
        % do input checks
        if any(cellfun(@(u) ~all(cellfun(@ismatrix, u)), U))
            error('ll1:U0', ['All U0{r}{n} should be matrices for n=1,2,4, and ' ...
                'a column vector for n=3.']);
        end
        if any(cellfun(@numel, U) ~= 4)
            error('ll1:U0', 'length(U0{r}) should be 4 for all r');
        end
        U = U(:).';
        U = cellfun(@(u) u(:).', U, 'UniformOutput', false);
        if any(cellfun(@(u) ~all(cellfun('size', u(1:3), 1)==size_tens),U))
            error('ll1:U0', ['size(U0{r}{n},1) should be size(T,n) for n=1:3 ' ...
                'and for all r']);
        end
        if length(U) ~= length(L)
            error('ll1:U0', 'length(U0) should be R=length(L)');
        end
        L1 = cellfun(@(u) size(u{1},2), U);
        L2 = cellfun(@(u) size(u{2},2), U);
        L3 = cellfun(@(u) size(u{3},2), U);
        if any(L(:) ~= L1(:)) || any(L(:) ~= L2(:))
            error('ll1:U0', 'size(U0{r}{n},2) should be L(r) for all r and n=1,2.');
        end
        if any(L3~=1)
            error('ll1:U0', 'size(U0{r}{3},2) should be 1 for all r.');
        end
        if any(cellfun(@(u) size(u{4},1), U)~=L) || ...
                any(cellfun(@(u) size(u{4},2), U)~=L)
            error('ll1:U0', ['U0{r}{4} should be a square L(r)xL(r) matrix for ' ...
                'all r.']);
        end
        for r = 1:length(U),
            U{r}{1} = U{r}{1}*U{r}{4};
            U{r}{4} = eye(size(U{r}{4}));
        end
    end
end

if strcmpi(options.OutputFormat, 'auto')
    if exist('U', 'var')
        options.OutputFormat = inputFormat;
    else
        options.OutputFormat = 'btd';
    end
end

if ~any(strcmpi(options.OutputFormat, {'cpd', 'btd'}))
    error('ll1:OutputFormat', 'Unknown output format %s', ...
        options.OutputFormat);
end

if ischar(options.Refinement) && strcmpi(options.Refinement, 'auto');
    options.Refinement = options.Algorithm;
end
if ~options.Display, print = @(varargin)true; else print = @fprintf; end
if ischar(options.ExploitStructure) && ~strcmp(options.ExploitStructure,'auto')
    error('ll1:exploitStructure','ExploitStructure option not supported!');
end

% STEP 0: Expand small structured tensors
expandlimitexceeded = prod(getsize(T)) > options.ExpandLimit;
output.Preprocessing.Name = 'none';
if isstructured && ~expandlimitexceeded
    output.Preprocessing.Name = 'expanded';
end

if isnumeric(options.Compression)
    if options.Compression == 0, options.Compression = false;
    elseif options.Compression == 1, options.Compression = true;
    else error('ll1:compression','Compression option not supported!');
    end
end

% STEP 1: compress the tensor if it was requested, and no initial factor
% matrices were provided by the user.
print('Step 1: Compression ');
if ~exist('U','var') && sum(L) > 1
    size_core = min(size_tens,[sum(L) sum(L), R]);
    ratio = prod(size_core)/prod(size_tens);
    if (ischar(options.Compression) && strcmpi(options.Compression,'auto'))|| ...
            (islogical(options.Compression) && options.Compression)
        if ratio < 0.5
            if any(strcmpi(type, {'full', 'sparse'}))
                options.Compression = @mlsvd_rsi;
                print('is mlsvd_rsi (compression ratio %.6g < 0.5)', ...
                    ratio);
            elseif ~strcmpi(type, 'incomplete') && prod(getsize(T)) <= ...
                    options.ExpandLimit
                options.Compression = @mlsvd_rsi;
                print('is mlsvd_rsi (compression ratio %.6g < 0.5)', ...
                    ratio);
            else
                options.Compression = @lmlra_aca;
                print('is lmlra_aca (compression ratio %.6g < 0.5)', ...
                    ratio);
            end
        else
            options.Compression = false;
            print('skipped (compression ratio %.6g >= 0.5).\n',ratio);
        end
    elseif xsfunc(options.Compression)
        print('is %s (requested by user)... ', func2str(options.Compression));
    else
        if options.Compression
            print('is mlsvd (requested by user).\n');
        else
            print('skipped (requested by user).\n');
        end
    end
elseif sum(L) == 1
    options.Compression = false;
    print('skipped (sum(L) = 1).\n');
else
    options.Compression = false;
    print('skipped (initialization U0 supplied).\n');
end

if xsfunc(options.Compression)
    if isstructured
        print(':\n');
        if ~expandlimitexceeded
            print('| Step 1a: Structured tensor expanded to full tensor (prod(getsize(T)) <= %.0g).\n', options.ExpandLimit);
            print(['| Step 1b: ' func2str(options.Compression) ' on full tensor... ']);
        else
            print('| Step 1a: Structured tensor expansion skipped (prod(getsize(T)) > %.0g).\n', options.ExpandLimit);
            print(['| Step 1b: ' func2str(options.Compression) ' on structured tensor... ']);
        end
    else
        print('... ');
    end
    
    Tunc = T;
    if ~isfield(options.CompressionOptions, 'FillCore')
        options.CompressionOptions.FillCore = true;
    end
    if isstructured && prod(getsize(T)) <= options.ExpandLimit
        T = ful(T);
    end
    [V,T] = options.Compression(T,size_core,options.CompressionOptions);
    output.Compression.Name = func2str(options.Compression);
    if sum(size_tens > 1) >= 3
        output.Compression.ratio = ratio;
        output.Compression.relerr = froblmlrares(Tunc,V,T)/frob(Tunc);
        print('relative error = %.6g.\n',output.Compression.relerr);
    else
        output.Compression.failed = true;
        options.Compression = false;
        T = Tunc; clear V;
        print('failed.\n');
    end
else
    if islogical(options.Compression) && options.Compression, print('\n'); end
    output.Compression.Name = 'none';
end

% STEP 2: initialize the factor matrices unless they were all provided by
% the user.
print('Step 2: Initialization ');
if ~exist('U','var'), U = cell(1,R); end
Uempty = cellfun(@isempty,U);
if any(Uempty)
    if ischar(options.Initialization) && ...
            strcmpi(options.Initialization,'auto')
        sz = getsize(T);
        sz = [sz ones(1, 3-length(sz))];
        if sz(1) >= sum(L) && sz(2) >= sum(L) && ...
                N == 3 && (sz(3) >= 2 || length(L) == 1)
            if isnumeric(T)
                options.Initialization = @ll1_gevd;
                print('is ll1_gevd (sum(sum(L) <= size(T)) >= 2)');
            elseif ~strcmpi(type, 'incomplete') && prod(sz) <= ...
                    options.ExpandLimit
                if ~exist('V', 'var')
                    Torig = T;
                    T = ful(T);
                end
                options.Initialization = @ll1_gevd;
                print('is ll1_gevd (sum(sum(L) <= size(T)) >= 2)');
            else
                options.Initialization = @ll1_rnd;
                print('is ll1_rnd (default)');
            end
        else
            options.Initialization = @ll1_rnd;
            print('is ll1_rnd (default)');
        end
    elseif isfunc(options.Initialization)
        print('is %s (requested by user)', ...
            func2str(options.Initialization));
    end
    if ~xsfunc(options.Initialization)
        error('ll1:Initialization','Not a valid initialization.');
    end
    
    % Make sure a complex solution is computed if requested
    if ~ischar(options.Complex) && options.Complex && isreal(ful(T,1));
        if strcmpi(func2str(options.Initialization), 'll1_rnd')
            if ~isfield(options.InitializationOptions, 'Imag')
                if isfield(options.InitializationOptions, 'Real')
                    options.InitializationOptions.Imag = ...
                        options.InitializationOptions.Real;
                else
                    options.InitializationOptions.Imag = @randn;
                end
            end
        elseif strcmpi(func2str(options.Initialization), 'll1_gevd')
            if ~isfield(options.InitializationOptions, 'isReal')
                options.InitializationOptions.isReal = false;
            end
        end
    end
else
    options.Initialization = false;
    print('is manual... ');
end
if isfunc(options.Initialization)
    if islogical(options.Compression) && ~options.Compression && isstructured
        print(':\n');
        if ~expandlimitexceeded
            print('| Step 2a: Structured tensor expanded to full tensor (prod(getsize(T)) <= %.0g).\n', options.ExpandLimit);
            print('| Step 2b: %s on full tensor... ',func2str(options.Initialization));
        else
            print('| Step 2a: Structured tensor expansion skipped (prod(getsize(T)) > %.0g).\n', options.ExpandLimit);
            print('| Step 2b: %s on structured tensor... ',func2str(options.Initialization));
        end
    else
        print('... ');
    end
    
    % Set the order in which the factor matrices should be updated.
    if sum(Uempty) > 0 && sum(Uempty) < N
        perm = [find(Uempty == 1) find(Uempty == 0)];
        if ~isfield(options,'AlgorithmOptions')
            options.AlgorithmOptions = struct;
        end
        if ~isfield(options.AlgorithmOptions,'Order')
            options.AlgorithmOptions.Order = perm;
        end
        if ~isfield(options.AlgorithmOptions, 'Display') && ...
                options.Display > 1
            options.AlgorithmOptions.Display = options.Display;
        end
    end
    % Generate initial factor matrices.
    [U0,output.Initialization] = options.Initialization(T,L,...
        options.InitializationOptions);
    % Fill empty factor matrices.
    for n = 1:sum(Uempty), U{n} = U0{n}; end
    output.Initialization.Name = func2str(options.Initialization);
    
    if exist('Torig', 'var')
        % reverse auto expansion of tensor
        T = Torig;
    end
else
    output.Initialization.Name = 'manual';
end
output.Initialization.relerr = frobll1res(T,U,L)/frob(T);
print('relative error = %.6g.\n',output.Initialization.relerr);

% STEP 3: run the selected LL1 algorithm.
if ~xsfunc(options.Compression) || (islogical(options.Refinement) && ...
        ~options.Refinement), unstring = 'un'; else unstring = ''; end

print(['Step 3: Algorithm is %s on the ' unstring 'compressed tensor'],func2str(options.Algorithm));

% STEP 3a: detect any structure in the tensor
structureDetected = false;

if xsfunc(options.Compression) && ischar(options.ExploitStructure) && strcmp(options.ExploitStructure,'auto')
    % print('| Step 3a: Detection of structure skipped (by default after compression).\n');
    print('... ');
else
    print(':\n');
    if ~options.ExploitStructure
        print('| Step 3a: Detection of structure skipped (requested by user).\n');
    else
        supportedStructures = {'hankel'};
        supportedStructuresPrint = {'Hankel'};
        
        if strcmp(getstructure(T),'full')
            print('| Step 3a: Detecting structure... ');
            for s = 1:numel(supportedStructures)
                [structure,representation] = detectstructure(T,'structures',supportedStructures{s});
                if ~isempty(structure)
                    print(['detected ' supportedStructuresPrint{s} ' structure.\n']);
                    structureDetected = true;
                    Tfull = T;
                    T = representation;
                    break
                end
            end
            if ~structureDetected, print('no structure detected.\n');
            else print('| Step 3b: Converted tensor to efficient representation.\n');
            end
        else
            print(['| Step 3a: Detection of structure skipped (' type ' structure already known).\n']);
        end
    end
end

% Step 3b: run the selected LL1 algorithm.
if xsfunc(options.Algorithm)
    if structureDetected
        % Structure detected
        print('| Step 3c: %s on efficient representation... ',func2str(options.Algorithm));
    elseif xsfunc(options.Compression) && ischar(options.ExploitStructure) && strcmp(options.ExploitStructure,'auto')
        % No detection of structure
    else
        % No structure detected
        print('| Step 3b: %s... ',func2str(options.Algorithm));
    end
    if ~isfield(options,'AlgorithmOptions')
        options.AlgorithmOptions = struct;
    end
    if ~isfield(options.AlgorithmOptions, 'TolX') && ~isnan(options.TolX)
        options.AlgorithmOptions.TolX = options.TolX;
    end
    if ~isfield(options.AlgorithmOptions, 'TolFun') && ~isnan(options.TolFun)
        options.AlgorithmOptions.TolFun = options.TolFun;
    end
    if ~isfield(options.AlgorithmOptions, 'MaxIter') && ~isnan(options.MaxIter)
        options.AlgorithmOptions.MaxIter = options.MaxIter;
    end
    if ~isfield(options.AlgorithmOptions, 'CGMaxIter') && ~isnan(options.CGMaxIter)
        options.AlgorithmOptions.CGMaxIter = options.CGMaxIter;
    end
    if ~isfield(options.AlgorithmOptions, 'Display') && options.Display > 1
        options.AlgorithmOptions.Display = options.Display;
    end
    
    warning('off','ll1_core:accuracy')
    UA0 = U;
    [U,output.Algorithm] = options.Algorithm(T,UA0,...
        options.AlgorithmOptions);
    output.Algorithm.Name = func2str(options.Algorithm);
else
    error('ll1:Algorithm','Not a valid algorithm.');
end
if isfield(output.Algorithm,'iterations')
    print('iterations = %i, ',output.Algorithm.iterations);
end
output.Algorithm.relerr = frobll1res(T,U,L)/frob(T);
print('relative error = %.6g.\n',output.Algorithm.relerr);

% STEP 3c: perform algorithm on full tensor
if structureDetected && output.Algorithm.info == 4
    print('| Step 3d: %s on the full tensor to improve accuracy... ',func2str(options.Algorithm));
    UB0 = U;
    [U,output.AlgorithmFull] = options.Algorithm(Tfull,UB0,...
        options.AlgorithmOptions);
    if isfield(output.AlgorithmFull,'iterations')
        print('iterations = %i, ',output.AlgorithmFull.iterations);
    end
    output.AlgorithmFull.relerr = frobll1res(Tfull,U)/frob(Tfull);
    print('relative error = %.6g.\n',output.AlgorithmFull.relerr);
    
    output.Algorithm = mergeOutputs(output.Algorithm,output.AlgorithmFull,UA0,UB0);
    output = rmfield(output,'AlgorithmFull');
end

% STEP 3d: expand to the original space.
if exist('V','var') && ~isempty(V)
    if all(cellfun(@iscell, U)) % BTD formatted
        for r = 1:R
            U{r}(1:N) = cellfun(@(u,v)v*u,U{r}(1:N),V,'UniformOutput',false);
        end
    else % CPD format
        U = cellfun(@(u,v) v*u, U, V, 'UniformOutput', false);
    end
    T = Tunc;
end

% STEP 4

refinement_skipped_nodecompression = ~xsfunc(options.Compression);
refinement_skipped_userrequested = islogical(options.Refinement) && ~options.Refinement;
continueRefinementThreshold = 1e-8;
refinement_skipped_thresholdstructuredreached = isstructured && (frob(ll1res(T,U))/frob(T) < continueRefinementThreshold);

% STEP 4a: detect any structure in the tensor
structureDetected = false;
if ~refinement_skipped_nodecompression && ~refinement_skipped_userrequested && ...
        ~refinement_skipped_thresholdstructuredreached
    print('Step 4: Refinement is %s on the uncompressed tensor:\n',func2str(options.Refinement));
    if options.ExploitStructure
        continueRefinementThreshold = 1e-8;
        if ~strcmp(type,'full')
            print(['| Step 4a: Detection of structure skipped (' type ' structure already known).\n']);
        else
            print('| Step 4a: Detecting structure... ');
            supportedStructures = {'hankel'};
            supportedStructuresPrint = {'Hankel'};
            for s = 1:numel(supportedStructures)
                [structure,representation] = detectstructure(Tunc,'structures',supportedStructures{s});
                if ~isempty(structure)
                    print(['detected ' supportedStructuresPrint{s} ' structure.\n']);
                    structureDetected = true;
                    Tuncfull = Tunc;
                    Tunc = representation;
                    break
                end
            end
            if ~structureDetected, print('no structure detected.\n');
            else print('| Step 4b: Converted tensor to efficient representation.\n');
            end
        end
    else
        %         if strcmp(unstring,'un'), unstring = ''; else unstring = 'un'; end
        print('| Step 4a: Detection of structure skipped (requested by user).\n');
    end
end

% STEP 4b: iteratively refine the solution if the tensor was compressed.
if ~xsfunc(options.Compression)
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (no decompression).\n');
elseif islogical(options.Refinement) && ~options.Refinement
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (requested by user).\n');
elseif refinement_skipped_thresholdstructuredreached
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (relative error <= %0.0g).\n',continueRefinementThreshold);
    if isfield(options.AlgorithmOptions,'TolFun') && options.AlgorithmOptions.TolFun <= eps
        warning('ll1:accuracyrefinement2',['Maximal numerical accuracy for'...
            ' efficient representations of structured tensors reached. The'...
            ' result may be improved using ll1 on ful(T) instead of on T.']);
    end
elseif xsfunc(options.Refinement)
    
    if structureDetected
        print('| Step 4c: %s on efficient representation... ',func2str(options.Refinement));
    else
        print('| Step 4b: %s... ',func2str(options.Refinement));
    end
    
    % Copy parameters from options.AlgorithmOptions which are not defined
    % in options.RefinementOptions.
    fn = fieldnames(options.AlgorithmOptions);
    for f = 1:length(fn)
        if ~isfield(options.RefinementOptions,fn{f})
            options.RefinementOptions.(fn{f}) = ...
                options.AlgorithmOptions.(fn{f});
        end
    end
    UA0 = U;
    [U,output.Refinement] = options.Refinement(Tunc,UA0, ...
        options.RefinementOptions);
    output.Refinement.Name = func2str(options.Refinement);
    if isfield(output.Refinement,'iterations')
        print('iterations = %i, ',output.Refinement.iterations);
    end
    output.Refinement.relerr = frobll1res(T,U,L)/frob(T);
    print('relative error = %.6g.\n',output.Refinement.relerr);
else
    error('ll1:Refinement','Not a valid refinement algorithm.');
end

% STEP 4c: perform algorithm on full tensor
if xsfunc(options.Refinement) && structureDetected && output.Refinement.info == 4
    print('| Step 4d: %s on the original full tensor to improve accuracy... ',func2str(options.Refinement));
    UB0 = U;
    [U,output.RefinementFull] = options.Refinement(Tuncfull,UB0,...
        options.RefinementOptions);
    if isfield(output.RefinementFull,'iterations')
        print('iterations = %i, ',output.RefinementFull.iterations);
    end
    output.RefinementFull.relerr = frobll1res(T,U)/frob(T);
    print('relative error = %.6g.\n',output.RefinementFull.relerr);
    
    output.Refinement = mergeOutputs(output.Refinement,output.RefinementFull,UA0,UB0);
    output = rmfield(output,'RefinementFull');
end

% Change output format if needed
if strcmpi(options.OutputFormat, 'cpd')
    tmp = U;
    U = cell(1, 3);
    for n = 1:3
        U{n} = cellfun(@(u) u{n}, tmp, 'UniformOutput', false);
        U{n} = cat(2, U{n}{:});
    end
end

% Format output.
fn = fieldnames(output);
for f = 1:length(fn)
    output.(fn{f}) = orderfields(output.(fn{f}));
end

    function output = mergeOutputs(outputA,outputB,UA0,UB0)
        % Merges outputs outputA and outputB, given initializations UA0 and
        % UB0 for first and second run, respectively.
        fieldsToKeep = {'Name','info','relerr'};
        outputAlgorithmFields = fieldnames(outputB);
        output = struct;
        for i = 1:numel(outputAlgorithmFields)
            field = outputAlgorithmFields{i};
            if isfield(outputA,field), valueA = outputA.(field);
            else valueA = []; end
            valueB = outputB.(field);
            switch field
                case 'fval'
                    % Merge values but omit final value from first run
                    output.(field) = [valueA(1:end-1) valueB];
                case 'iterations'
                    % Sum values from both runs
                    output.(field) = valueA + valueB;
                case 'relfval'
                    % Merge values but renormalize values from second run
                    output.(field) = cat(2,valueA,valueB*...
                        abs(outputB.fval(1))/...
                        abs(outputA.fval(1)));
                case 'relstep'
                    % Merge values but renormalize values from second run
                    UB0 = ll1convert(UB0);
                    UA0 = ll1convert(UA0);
                    output.(field) = cat(2,valueA,valueB*...
                        sqrt(sum(cellfun(@(x) frob(x,'squared'),UB0)))/...
                        sqrt(sum(cellfun(@(x) frob(x,'squared'),UA0))));
                otherwise
                    if ischar(valueB) || any(strcmp(field,fieldsToKeep)) || numel(valueB)==1
                        % Take value from second run
                        output.(field) = valueB;
                    else
                        % Merge values from both runs
                        output.(field) = cat(2,valueA,valueB);
                    end
            end
        end
    end
end
