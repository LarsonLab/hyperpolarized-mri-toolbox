function [U,output] = cpd(T,R,varargin)
%CPD Canonical polyadic decomposition.
%   [U,output] = CPD(T,R) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T in R rank-one terms.
%
%   [U,output] = CPD(T,U0) allows the user to provide initial factor
%   matrices U0, to be used by the main algorithm. By providing an
%   initialization, the compression step will be skipped. Factor matrices
%   equal to the empty matrix [] will be initialized using the chosen
%   initialization method.
%
%   The structure output contains the output of the selected algorithms:
%
%      output.Preprocessing  - The output of the preprocessing step
%      output.Compression    - The output of the compression step.
%      output.Initialization - The output of the initialization step.
%      output.Algorithm      - The output of the main algorithm.
%      output.Refinement     - The output of the the refinement step after
%                              decompression (if applicable).
%
%   CPD(T,R,options) and CPD(T,U0,options) allow the user to choose which
%   algorithms will be used in the different steps and also set parameters
%   for those algorithms:
%
%      Use options.Compression = [{'auto'}|@mlsvd_rsi,@lmlra_aca,@mlsvd] to
%      select whether or not the tensor is first compressed using the given
%      algorithm to a tensor of dimensions min(size(T),R*ones(1,N)). If one
%      or more initial factor matrices are provided compression is skipped.
%      On 'auto', compression is only executed if it is expected to lead to
%      a significant reduction in computational complexity. In the case the
%      tensor is full or sparse, or if the tensor is not incomplete and
%      prod(getsize(T)) <= options.ExpandLimit (see further), mlsvd_rsi is
%      the default compression algorithm for 'auto'; in all other cases the
%      default is lmlra_aca. The structure options.CompressionOptions will
%      be passed to the chosen algorithm.
%
%      Use options.Initialization = [{'auto'}|@cpd_rnd|@cpd_gevd] to
%      choose the initialization method. The structure
%      options.InitializationOptions will be passed to the chosen
%      initialization method. On 'auto', @cpd_gevd is used when possible,
%      otherwise @cpd_rnd.
%
%      Use options.Algorithm = [@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|...
%      @cpd3_sgsd|@cpd_rbs] to choose the main algorithm. The structure
%      options.AlgorithmOptions will be passed to the chosen algorithm.
%
%      Use options.Refinement = [@cpd_als|@cpd_minf|{@cpd_nls}|@cpd3_sd|...
%      @cpd3_sgsd|cpdnn_nlsb|cpd_rbs|false] to choose the algorithm used to
%      refine the solution after decompression, if applicable. By default,
%      options.Refinement is equal to options.Algorithm. The structure
%      options.RefinementOptions will be passed to the chosen refinement
%      method. Parameters not set in options.RefinementOptions will be
%      copied from options.AlgorithmOptions where possible.
%
%   Further options are:
%
%      options.Display = false   - Set to true to enable printing output
%                                  information to the command line. If
%                                  options.Display > 1, it is passed on to
%                                  AlgorithmOptions and RefinementOptions,
%                                  unless these structs define Display.
%      options.ExpandLimit =     - Create full tensor if T is not incomplete
%       1e6                        and prod(getsize(T)) <=
%                                  options.ExpandLimit for compression
%                                  and/or initialization steps
%      options.Complex =         - If true, compute a complex solution even if
%      ['auto', true]              the tensor is real. If auto, a complex
%                                  solution is only computed if T is
%                                  complex,  or if options.Initialization
%                                  results in a complex initialization,
%                                  e.g., if options.InitializationOptions
%                                  .Imag = @randn. This option only works
%                                  for options.Initialization = cpd_rnd or
%                                  cpd_gevd.
%      options.ExploitStructure  - Detect structure in the tensor and
%       = [false,true,{'auto'}]    exploit this for speed purposes. If
%                                  true, this is applied for the
%                                  decompositions of both the compressed
%                                  and uncompressed tensor. If auto, this
%                                  is only applied for the decomposition of
%                                  the uncompressed tensor.
%
%   The following options are passed to AlgorithmOptions and
%   RefinementOptions unless these structs define these options:
%      options.TolX            - Tolerance for step length
%      options.TolFun          - Tolerance for function value
%      options.MaxIter         - Maximum number of iterations
%      options.CGMaxIter       - Maximum number of CG iterations for inexact
%                                NLS type algorithms.
%      options.UseCPDI         - If true, the CPDI routines are used for
%                                incomplete tensors
%
%   See also rankest, cpdgen, cpderr.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/03/01   OD      Added structure detection
% - 2016/01/28   NV      Added structured tensors

% Check the tensor T.
type = getstructure(T);
isstructured = ~any(strcmp(type, {'full', 'sparse', 'incomplete'}));

if ~isstructured
    T = fmt(T);
    type = getstructure(T);
end
size_tens = getsize(T);
N = length(size_tens);

if iscell(R)
    U = R;
    R = max(cellfun('size',U,2));
    N = length(U);
    
    sz = cellfun('size', U, 1);
    if length(sz) < length(size_tens)
        error('cpd:U0', ['N = length(U0) should be >= getorder(T). If U0{n} should ' ...
            'be automatically initialized for n > N, set ' ...
            'U0{n} = [].']);
    end
    if all(cellfun(@isempty, U))
        error('cpd:U0', ['Could not determine R as all U0{n} are empty. Use ' ...
            'cpd(T,R) instead or U0{n} should be a matrix for at ' ...
            'least one n.']);
    end
    if any(sz(1:length(size_tens)) ~= size_tens & sz(1:length(size_tens))>0)
        error('cpd:U0', ['size(U0{n},1) should be size(T,n) or U0{n} should ' ...
            'be empty for all n']);
    end
    if any(sz(length(size_tens)+1:end)~=1)
        error('cpd:U0', 'size(U0{n},1) should be 1 for all n > getorder(T)');
    end
    if any(cellfun('size',U,2) ~= R & cellfun('size',U,2)~=0)
        error('cpd:U0', 'size(U0{n},2) should be R for all n');
    end
end

if isstructured
    try 
        frob(T);
    catch e
        if strcmpi(e.identifier, 'frob:notImplemented');
            error('cpd:notImplemented', ...
                  ['cpd does not support the structured tensor type %s, yet. Use ' ...
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
p.addOptional('Algorithm', @cpd_nls); %@cpd_nls,@cpd_als,@cpd_minf,@cpd3_sd,@cpd3_sgsd
p.addOptional('AlgorithmOptions', struct());
p.addOptional('Refinement', 'auto');
p.addOptional('RefinementOptions', struct());
p.addOptional('Display', false);
p.addOptional('Complex', 'auto');
p.addOptional('ExpandLimit', 1e6);
p.addOptional('TolX', nan);
p.addOptional('TolFun', nan);
p.addOptional('MaxIter', nan);
p.addOptional('CGMaxIter', nan);
p.addOptional('UseCPDI', nan);
p.addOptional('ExploitStructure','auto');
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

if ischar(options.Refinement) && strcmpi(options.Refinement, 'auto');
    options.Refinement = options.Algorithm;
end
if ~options.Display, print = @(varargin)true; else print = @fprintf; end
if ischar(options.ExploitStructure) && ~strcmp(options.ExploitStructure,'auto')
    error('cpd:exploitStructure','ExploitStructure option not supported!');
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
    else error('cpd:compression','Compression option not supported!');
    end
end

% STEP 1: compress the tensor if it was requested, and no initial factor
% matrices were provided by the user.
print('Step 1: Compression ');
if ~exist('U','var') && R > 1
    size_core = min(size_tens,R*ones(1,N));
    ratio = prod(size_core)/prod(size_tens);
    
    % If compression is set to auto ('auto' or true).
    if (ischar(options.Compression) && strcmpi(options.Compression,'auto')) || ...
            (islogical(options.Compression) && options.Compression)
        if ratio < 0.5
            if any(strcmp(type, {'full', 'sparse'}))
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
    % Instead, if a custom function handle is given for compression
    elseif xsfunc(options.Compression)
        print('is %s (requested by user)', func2str(options.Compression));
    else
        % If compression is set to 'mlsvd'
        if options.Compression
            print('is mlsvd (requested by user)');
        else
            print('skipped (requested by user).\n');
        end
    end
elseif R == 1
    options.Compression = false;
    print('skipped (R = 1).\n');
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
    if sum(size(T) > 1) >= 3
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
if ~exist('U','var'), U = cell(1,N); end
Uempty = cellfun(@isempty,U);
if any(Uempty)
    if ischar(options.Initialization) && ...
            strcmpi(options.Initialization,'auto')
        sz = getsize(T);
        sz = [sz ones(1, N-length(sz))];
        if sum(R <= sz) >= 2 && N > 2
            if isnumeric(T)
                options.Initialization = @cpd_gevd;
                print('is cpd_gevd (sum(R <= size(T)) >= 2)');
            elseif ~strcmpi(type, 'incomplete') && prod(sz) <= ...
                    options.ExpandLimit
                if ~exist('V', 'var')
                    Torig = T;
                    T = ful(T);
                end
                options.Initialization = @cpd_gevd;
                print('is cpd_gevd (sum(R <= size(T)) >= 2)');
            else
                options.Initialization = @cpd_rnd;
                print('is cpd_rnd (default)');
            end
        else
            options.Initialization = @cpd_rnd;
            print('is cpd_rnd (default)');
        end
    elseif isfunc(options.Initialization)
        print('is %s (requested by user)', ...
            func2str(options.Initialization));
    end
    if ~xsfunc(options.Initialization)
        error('cpd:Initialization','Not a valid initialization.');
    end
    
    % Make sure a complex solution is computed if requested
    if ~ischar(options.Complex) && options.Complex && isreal(ful(T,1))
        if strcmpi(func2str(options.Initialization), 'cpd_rnd')
            if ~isfield(options.InitializationOptions, 'Imag')
                if isfield(options.InitializationOptions, 'Real')
                    options.InitializationOptions.Imag = ...
                        options.InitializationOptions.Real;
                else
                    options.InitializationOptions.Imag = @randn;
                end
            end
        elseif strcmpi(func2str(options.Initialization), 'cpd_gevd')
            if ~isfield(options.InitializationOptions, 'IsReal')
                options.InitializationOptions.IsReal = false;
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
    % Add orthogonality option for cpd_rnd
    if strcmp(func2str(options.Initialization), 'cpd_rnd')
        if ~isfield(options.InitializationOptions, 'Orth') 
            options.InitializationOptions.Orth = true;
        end
    end
    
    % Generate initial factor matrices.
    [U0,output.Initialization] = options.Initialization(T,R,...
        options.InitializationOptions);
    
    % Fill empty factor matrices.
    if sum(Uempty) > 0 && sum(Uempty) < length(Uempty)
        [~,P,D] = cpderr(U(~Uempty), U0(~Uempty));
        U0(Uempty) = cellfun(@(u) u*P, U0(Uempty), 'UniformOutput', false);
        D = cellfun(@diag, D,'UniformOutput', false);
        D = prod(cat(2,D{:}),2);
        U0{find(Uempty, 1)} = U0{find(Uempty,1)}*diag(1./D);
    end
    U(Uempty) = U0(Uempty);
    output.Initialization.Name = func2str(options.Initialization);
    
    if exist('Torig', 'var')
        % reverse auto expansion of tensor
        T = Torig;
    end
else
    output.Initialization.Name = 'manual';
end
output.Initialization.relerr = frobcpdres(T,U)/frob(T);
print('relative error = %.6g.\n',output.Initialization.relerr);

% STEP 3

if ~xsfunc(options.Compression) || (islogical(options.Refinement) && ...
        ~options.Refinement)
    unstring = 'un';
else
    unstring = '';
end

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


% STEP 3b: run the selected CPD algorithm.
if xsfunc(options.Algorithm)
    if structureDetected
        % Structure detected
        print('| Step 3c: %s on the efficient representation... ',func2str(options.Algorithm));
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
    if ~isfield(options.AlgorithmOptions, 'UseCPDI') && ~isnan(options.UseCPDI)
        options.AlgorithmOptions.UseCPDI = options.UseCPDI;
    end
    if ~isfield(options.AlgorithmOptions, 'Display') && options.Display > 1
        options.AlgorithmOptions.Display = options.Display;
    end
    useincompletecore = false;
    if (~isstruct(T) || (isfield(T,'incomplete') && ~T.incomplete)) && ...
            isfield(options.AlgorithmOptions, 'UseCPDI')
        options.AlgorithmOptions.UseCPDI = false;
        useincompletecore = true;
    end
    
    warning('off','cpd_core:accuracy');
    UA0 = U;
    [U,output.Algorithm] = options.Algorithm(T,UA0,...
        options.AlgorithmOptions);
    output.Algorithm.Name = func2str(options.Algorithm);
    if useincompletecore
        options.AlgorithmOptions.UseCPDI = true;
    end
else
    error('cpd:Algorithm','Not a valid algorithm.');
end
if isfield(output.Algorithm,'iterations')
    print('iterations = %i, ',output.Algorithm.iterations);
end
output.Algorithm.relerr = frobcpdres(T,U)/frob(T);
print('relative error = %.6g.\n',output.Algorithm.relerr);

warning('on','cpd_core:accuracy');
if (~xsfunc(options.Compression) || (islogical(options.Refinement) && ...
        ~options.Refinement)) && output.Algorithm.info==4 && isstructured
    warning('cpd:accuracyalgorithm',['Maximal numerical accuracy for'...
        ' structured tensors reached. The result may be improved using' ...
        ' cpd on ful(T) instead of on T.']);
end
    

% STEP 3c: perform algorithm on full tensor
if structureDetected && output.Algorithm.info == 4
    print('| Step 3d: %s on the full tensor to improve accuracy... ',func2str(options.Algorithm));
    UB0 = U;
    [U,output.AlgorithmFull] = options.Algorithm(Tfull,UB0,...
        options.AlgorithmOptions);
    if isfield(output.AlgorithmFull,'iterations')
        print('iterations = %i, ',output.AlgorithmFull.iterations);
    end
    output.AlgorithmFull.relerr = frobcpdres(Tfull,U)/frob(Tfull);
    print('relative error = %.6g.\n',output.AlgorithmFull.relerr);
    
    output.Algorithm = mergeOutputs(output.Algorithm,output.AlgorithmFull,UA0,UB0);
    output = rmfield(output,'AlgorithmFull');
end

% STEP 3d: expand to the original space.
if exist('V','var') && ~isempty(V)
    U = cellfun(@(u,v)v*u,U,V,'UniformOutput',false);
    T = Tunc;
end

% STEP 4

refinement_skipped_nodecompression = ~xsfunc(options.Compression);
refinement_skipped_userrequested = islogical(options.Refinement) && ~options.Refinement;
% refinement_skipped_TolFunreached = isfield(options.RefinementOptions,'TolFun') && .5*(frob(cpdres(T,U))^2 < options.RefinementOptions.TolFun);
continueRefinementThreshold = 1e-8;
refinement_skipped_thresholdstructuredreached = isstructured && (frob(cpdres(T,U))/frob(T) < continueRefinementThreshold);

% STEP 4a: detect any structure in the tensor

structureDetected = false;
if ~refinement_skipped_nodecompression && ~refinement_skipped_userrequested && ...
        ~refinement_skipped_thresholdstructuredreached
    print('Step 4: Refinement is %s on the uncompressed tensor:\n',func2str(options.Refinement));
    if options.ExploitStructure
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
if refinement_skipped_nodecompression
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (no decompression).\n');
elseif refinement_skipped_userrequested
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (requested by user).\n');
% elseif refinement_skipped_TolFunreached
%     output.Refinement.Name = 'none';
%     print('Step 4: Refinement ');
%     print('skipped (relative error ).\n');
elseif refinement_skipped_thresholdstructuredreached
    output.Refinement.Name = 'none';
    print('Step 4: Refinement ');
    print('skipped (relative error of %.6g <= %0.0g).\n',frob(cpdres(T,U))/frob(T),...
        continueRefinementThreshold);
    if isfield(options.AlgorithmOptions,'TolFun') && options.AlgorithmOptions.TolFun <= eps
        warning('cpd:accuracyrefinement2',['Maximal numerical accuracy for'...
            ' efficient representations of structured tensors reached. The'...
            ' result may be improved using cpd on ful(T) instead of on T.']);
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
    if (~isstruct(T) || (isfield(T,'incomplete') && ~T.incomplete)) && ...
            isfield(options.RefinementOptions, 'UseCPDI')
        % disable cpdi if not an incomplete tensor
        options.RefinementOptions.UseCPDI = false;
    end
    
    warning('off','cpd_core:accuracy');
    UA0 = U;
    [U,output.Refinement] = options.Refinement(Tunc,UA0, ...
        options.RefinementOptions);
    output.Refinement.Name = func2str(options.Refinement);
    if isfield(output.Refinement,'iterations')
        print('iterations = %i, ',output.Refinement.iterations);
    end
    output.Refinement.relerr = frobcpdres(T,U)/frob(T);
    print('relative error = %.6g.\n',output.Refinement.relerr);
    warning('on','cpd_core:accuracy');
    if output.Refinement.info==4 && isstructured
        warning('cpd:accuracyrefinement',['Maximal numerical accuracy for'...
            ' structured tensors reached. The result may be improved using' ...
            ' cpd on ful(T) instead of on T.']);
    end
else
    error('cpd:Refinement','Not a valid refinement algorithm.');
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
    output.RefinementFull.relerr = frobcpdres(T,U)/frob(T);
    print('relative error = %.6g.\n',output.RefinementFull.relerr);
    
    output.Refinement = mergeOutputs(output.Refinement,output.RefinementFull,UA0,UB0);
    output = rmfield(output,'RefinementFull');
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

