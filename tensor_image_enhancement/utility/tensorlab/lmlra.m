function [U,S,output] = lmlra(T,size_core,varargin)
%LMLRA Low multilinear rank approximation.
%   [U,S,output] = lmlra(T,size_core) computes the factor matrices U{1},
%   ..., U{N} and core tensor S of dimensions size_core belonging to a low
%   multilinear rank approximation of the N-th order tensor T.
%
%   [U,S,output] = lmlra(T,U0,S0) allows the user to provide initial factor
%   matrices U0 and core tensor S0, to be used by the main algorithm. Factor
%   matrices equal to the empty matrix [] will be initialized using the chosen
%   initialization method. If exactly one factor matrix n is empty and S0 is
%   given, the empty matrix is initialized using a least squares solution, if
%   one exists with size(S0,n) columns. If S0 is not given or empty, S0 is
%   computed as tmprod(T, V, 1:N) in which V{n} = pinv(U{n}) and U{n} = U0{n}
%   if given, or determined by the initalization method. 
%
%   The structure output contains the output of the selected algorithms:
%
%      output.Preprocessing  - The output of the preprocessing step
%      output.Initialization - The output of the initialization step.
%      output.Algorithm      - The output of the main algorithm.
%
%   lmlra(T,size_core,options) and lmlra(T,U0,S0,options) allow the user to
%   choose which algorithms will be used in the different steps and also
%   set parameters for those algorithms:
%
%      Use options.Initialization = [{'auto'}|@mlsvd_rsi|@lmlra_aca|@lmlra_rnd|
%      @mlsvd] to choose the initialization method. By default the randomized
%      mlsvd_rsi method is used if the tensor is full or spare, or if the tensor
%      is not incomplete and prod(getsize(T)) <= options.ExpandLimit (see
%      further). For all other tensors an adaptive cross approximation algorithm
%      is used. The struct options.InitializationOptions will be passed to the
%      chosen initialization method as third argument.
%
%      Use options.Algorithm = [@lmlra_hooi|@lmlra_minf|{@lmlra_nls}| ...
%      lmlra3_dgn|lmlra3_rtr] to choose the main algorithm. The structure
%      options.AlgorithmOptions will be passed to the chosen algorithm.
%
%   Further options are:
%
%      options.Normalize = true  - Normalize the output of the LMLRA
%                                  (orthogonal factors U, ordered,
%                                  all-orthogonal core tensor S).   
%      options.Display = false   - Set to true to enable printing output
%                                  information to the command line. If
%                                  options.Display > 1, it is passed on to
%                                  AlgorithmOptions and RefinementOptions,
%                                  unless these structs define Display.
%      options.ExpandLimit =     - Create full tensor if T is not incomplete
%      1e6                         and prod(getsize(T)) <= options.ExpandLimit
%                                  for compression and/or initialization
%                                  steps
%
%   The following options are passed on to AlgorithmOptions unless these
%   structs define these options:
%      options.TolX            - Tolerance for step length
%      options.TolFun          - Tolerance for function value
%      options.MaxIter         - Maximum number of iterations
%      options.CGMaxIter       - Maximum number of CG iterations for inexact
%                                NLS type algorithms.
%
%   See also mlrankest, lmlragen, lmlraerr.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Gather some input information.
type = getstructure(T);
isstructured = ~any(strcmp(type, {'full', 'sparse', 'incomplete'}));
if ~isstructured
    T = fmt(T);
    type = getstructure(T);
end
size_tens = getsize(T);
N = length(size_tens);

if nargin < 3, S0 = [];  end
if nargin >= 3 && (isnumeric(varargin{1}) || isempty(varargin{1}))
    % S0 is given
    S = varargin{1};
    S0 = S;
    varargin = varargin(2:end);
end
if iscell(size_core)
    % U0 is given
    U = size_core(:).';
    size_core = cellfun('size', U, 2);
    if ~exist('S0', 'var')
        S0 = [];
    end
end    
if isempty(size_core)
    U = cell(1, N);
end
if isempty(size_core) || any(size_core == 0)
    if exist('S', 'var') && ~isempty(S)
        size_core = size(S);
    else 
        error('lmlra:size_core', ['Could not determine size_core. If size(U0{n},2) ' ...
                            '== 0 for some n, S0 should be given and should ' ...
                            'be non-empty.'])
    end    
end
size_core = size_core(:).';

% Check inputs
if exist('U', 'var')
    sz = cellfun('size', U, 1);
    if length(sz) < length(size_tens)
        error('lmlra:U0', 'length(U) should be >= getorder(T)');
    end
    if exist('S', 'var') && ~isempty(S)
        size_S = size(S0);
        if length(size_S) > length(size_core)
            error('lmlra:U0', 'length(U0) should be >= ndims(S0)');
        end
        if any(size_core(1:length(size_S)) ~= size_S)
            error('lmlra:U0', 'size(U0{n},2) should be size(S0,n) for all n');
        end
        if any(size_core(length(size_S)+1:end)~=1) 
            error('lmlra:U0', 'size(U0{n},2) should be 1 for all n > ndims(S0)');
        end
    end 
    if any(sz(1:length(size_tens)) ~= size_tens & sz(1:length(size_tens)) > 0)
        error('lmlra:U0', 'size(U0{n},1) should be size(T,n) for all n');
    end
    if any(sz(length(size_tens)+1:end)~=1) 
        error('lmlra:U0', 'size(U0{n},1) should be 1 for all n > getorder(T)');
    end
    size_core(length(size_core)+1:length(U)) = 1;
end 

try 
    frobT = frob(T);
catch e
    if strcmpi(e.identifier, 'frob:notImplemented');
        error('lmlra:notImplemented', ...
              ['lmlra does not support the structured tensor type %s, yet. Use ' ...
               'ful(T) instead.'], type);
    end
end

% Check the options structure.
p = inputParser();
p.addOptional('Normalize', true);
p.addOptional('Initialization', 'auto');
p.addOptional('InitializationOptions', struct());
p.addOptional('Algorithm', @lmlra_nls);
p.addOptional('AlgorithmOptions', struct());
p.addOptional('ExpandLimit', 1e6);
p.addOptional('TolX', nan);
p.addOptional('TolFun', nan);
p.addOptional('MaxIter', nan);
p.addOptional('CGMaxIter', nan);
p.addOptional('Display', false);
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if length(size_core) == 1, options.Initialization = @mlsvd; end
if ~options.Display, print = @(varargin)true; else print = @fprintf; end

% Set default preprocessing to none (to be overwritten later)
output.Preprocessing.Name = 'none';

% STEP 1: initialize the factor matrices unless they were all provided by
% the user.
if ischar(options.Initialization) 
    if strcmp(options.Initialization, 'auto')
        if any(strcmp(type, {'full', 'sparse'}))
            options.Initialization = @mlsvd_rsi;
        elseif ~strcmpi(type, 'incomplete') && prod(getsize(T)) <= ...
                options.ExpandLimit
            options.Initialization = @mlsvd_rsi;     
        else 
            options.Initialization = @lmlra_aca;
        end
    else 
        error('lmlra:Initialization', ['The initialization function should be ' ...
                            'a function handle, or the string ''auto''.'])
    end
end

print('Step 1: Initialization ');
if ~exist('U','var'), U = cell(1,N); end
Uempty = cellfun(@isempty,U);
if any(Uempty) || isempty(S0)
    if isfunc(options.Initialization)
        print('is %s',func2str(options.Initialization));
    end
    if ~xsfunc(options.Initialization)
        error('lmlra:Initialization','Not a valid initialization.');
    end
else  
    options.Initialization = false;
    print('is manual... ');
end
if isfunc(options.Initialization)
    
    if isstructured && prod(getsize(T)) <= options.ExpandLimit
        output.Preprocessing.Name = 'expanded';
        print([':\n| Step 1a: Structured tensor expanded to full tensor ' ...
               '(prod(getsize(T) <= %.0g)\n'], options.ExpandLimit);
        Torig = T;
        T = ful(T);
        print('| Step 1b: %s on full tensor... ', ...
              func2str(options.Initialization));
    else 
        print('...');
    end
    
    projectS0 = false;
    if sum(Uempty) == 0 && isempty(S0) 
        % If only core missing, find core using given factors
        S0 = tmprod(T, cellfun(@pinv,U,'uni',0), 1:length(U));
    elseif sum(Uempty) == 1 && ~isempty(S0)
        % If one factor missing, solve least squares problem
        idx = find(~Uempty);
        tmp = tmprod(T, cellfun(@pinv,U(idx),'uni',0), idx);
        U0 = U;
        n = find(Uempty);
        U0{Uempty} = tens2mat(tmp,n)/tens2mat(S0,n);
    else 
        % Generate initial factor matrices.
        if sum(Uempty) ~= length(Uempty) && (~exist('S0','var') || isempty(S0))
            projectS0 = true; 
        end            
        [U0,S0,sv] = options.Initialization(T,size_core,...
                                            options.InitializationOptions);
        if isstruct(sv)
            output.Initialization = sv; 
        else
            output.Initialization.sv = sv;
        end    
    end 

    % Fill empty factor matrices.
    if sum(Uempty) > 0 
        U(Uempty) = U0(Uempty);
    end
    if projectS0 && sum(Uempty) > 0
        Utmp = cellfun(@pinv, U, 'UniformOutput', false);
        S = reshape(full(mtkronprod(T, Utmp, 0, 'H')), size_core);
    elseif ~exist('S', 'var') || isempty(S) 
        S = S0;
    end
    output.Initialization.Name = func2str(options.Initialization);
else 
    output.Initialization.Name = 'manual';
end
output.Initialization.relerr = froblmlrares(T,U,S)/frob(T);
print('relative error = %.6g.\n',output.Initialization.relerr);

if strcmpi(output.Preprocessing.Name, 'expanded')
    T = Torig;
end

% STEP 2: run the main LMLRA algorithm.
if xsfunc(options.Algorithm)
    print('Step 2: Algorithm is %s... ',func2str(options.Algorithm));
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
    [U,S,output.Algorithm] = options.Algorithm(T,U,S,...
                             options.AlgorithmOptions);
    output.Algorithm.Name = func2str(options.Algorithm);
else
    error('lmlra:Algorithm','Not a valid algorithm.');
end
if isfield(output.Algorithm,'iterations')
    print('iterations = %i, ',output.Algorithm.iterations);
end
output.Algorithm.relerr = froblmlrares(T,U,S)/frobT;
print('relative error = %.6g.\n',output.Algorithm.relerr);

% Format output.
fn = fieldnames(output);
for f = 1:length(fn)
    output.(fn{f}) = orderfields(output.(fn{f}));
end

% normalize results
if options.Normalize
    if any(cellfun(@(u) frob(u'*u - eye(size(u,2))), U) > 1e-14)
        [U,R] = cellfun(@(u) qr(u, 0), U, 'UniformOutput', false);
        S = lmlragen(R,S);
    end
    [tmp, S] = mlsvd(S);
    U(1:numel(tmp)) = cellfun(@(u,v) u*v, U(1:numel(tmp)), tmp, 'UniformOutput', false);
end 
