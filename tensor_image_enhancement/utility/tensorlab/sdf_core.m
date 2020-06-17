function [sol,output] = sdf_core(model,varargin)
%SDF_CORE Computational core for structured data fusion.
%   SDF_CORE should not be called directly. Use SDF_MINF or SDF_NLS instead. 
%
%   See also: sdf_minf, sdf_nls

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%
% Version History:
% - 2016/01/07   NV      Added regL0
% - 2016/12/14   NV      Added ll1 
% - 2016/12/12   NV      Added cpdi 
% - 2015/12/10   NV      Faster reference resolution, skipped constants,
%                        bypasses for untransformed factors.
        

% Check the options structure.
p = inputParser;
p.addOptional('OptimizationType', 'minf');
p.addOptional('MaxIter', 5000);
p.addOptional('CGMaxIter', 15);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-8);
p.addOptional('TolAbs', 0);
p.addOptional('TolLargeScale', 0.02);
p.addOptional('JHasFullRank', false);
p.addOptional('Display', 0);
p.addOptional('PC', true);
p.addOptional('CheckModel', true);
p.addOptional('Weights', []);
p.addOptional('RelWeights', []);
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

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
        error('sdf:incompatibleAlgorithm', ...
              ['The %s method is incompatible with the nls optimization ' ...
               'type.'], func2str(options.Algorithm));
    elseif ~ismember(func2str(options.Algorithm), nlsfun)
        warning('sdf:unknownAlgorithm', ['The %s is not known by SDF. Use ' ...
                            'at your own risk'], options.Algorithm);
    end
  case 'minf'
    if ~isfield(options, 'Algorithm')
        options.Algorithm = @minf_lbfgsdl;
    elseif ismember(func2str(options.Algorithm), nlsfun)
        error('sdf:incompatibleAlgorithm', ...
              ['The %s method is incompatible with the minf optimization ' ...
               'type.'], func2str(options.Algorithm));
    elseif ~ismember(func2str(options.Algorithm), minffun)
        warning('sdf:unknownAlgorithm', ['The %s is not known by SDF. Use ' ...
                            'at your own risk'], options.Algorithm);
    end
  otherwise
    error('sdf:unknownOptimizationType', ...
          ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
           'are supported now'], options.OptimizationType);
end
        
% Check consistency of model
try 
    if options.CheckModel 
        model = sdf_check(model, 'internal', true);
    end
catch e
    rethrow(e)
end

% Set functions for saving state, computing the objective function,
% gradient and fast matrix-vector products. New models need only implement
% these four functions.
for I = 1:length(model.factorizations)
    switch model.factorizations{I}.type
      case 'cpd'
        if model.factorizations{I}.isstructured 
            model.factorizations{I}.objfun = @objfun_cpd_structured;
            model.factorizations{I}.grad = @grad_cpd_structured;
        else 
            model.factorizations{I}.objfun = @objfun_cpd;
            model.factorizations{I}.grad = @grad_cpd;
        end 
        model.factorizations{I}.state = @state_cpd;
        model.factorizations{I}.JHJx = @JHJx_cpd;
      case 'cpdi'
        model.factorizations{I}.objfun = @objfun_cpd;
        model.factorizations{I}.grad = @grad_cpd;
        model.factorizations{I}.state = @state_cpd;
        model.factorizations{I}.JHJx = @JHJxi_cpd;  
      case 'll1'
        if model.factorizations{I}.isstructured 
            model.factorizations{I}.objfun = @objfun_ll1_structured;
            model.factorizations{I}.grad = @grad_ll1_structured;
        else 
            model.factorizations{I}.objfun = @objfun_ll1;
            model.factorizations{I}.grad = @grad_ll1;
        end 
        model.factorizations{I}.state = @state_ll1;
        model.factorizations{I}.JHJx = @JHJx_ll1;
      case {'btd', 'lmlra'}
        if model.factorizations{I}.isstructured 
            model.factorizations{I}.objfun = @objfun_btd_structured;
            model.factorizations{I}.grad = @grad_btd_structured;
        else
            model.factorizations{I}.objfun = @objfun_btd;
            model.factorizations{I}.grad = @grad_btd;
        end
        model.factorizations{I}.state = @state_btd;
        model.factorizations{I}.JHJx = @JHJx_btd;
      case 'regL2'
        model.factorizations{I}.state = @state_regL2;
        model.factorizations{I}.objfun = @objfun_regL2;
        model.factorizations{I}.grad = @grad_regL2;
        model.factorizations{I}.JHJx = @JHJx_regL2;
      case 'regL1'
        model.factorizations{I}.state = @state_regL1;
        model.factorizations{I}.objfun = @objfun_regL1;
        model.factorizations{I}.grad = @grad_regL1;
        model.factorizations{I}.JHJx = @JHJx_regL1;
      case 'regL0'
        model.factorizations{I}.state = @state_regL0;
        model.factorizations{I}.objfun = @objfun_regL0;
        model.factorizations{I}.grad = @grad_regL0;
        model.factorizations{I}.JHJx = @JHJx_regL0;
    end
end

% Fill in constants and dereference pointers from factors to variables.
cache.factors.isconst = model.isconstant; %TODO remove constant from cache

% Convert full data to incomplete/sparse format where appropriate.
cache.factorizations.isincomplete = false(1,length(model.factorizations));
cache.factorizations.issparse = false(1,length(model.factorizations));
for I = 1:length(model.factorizations)
    if ~isfield(model.factorizations{I},'data') || ...
       iscell(model.factorizations{I}.data)
        continue;
    end
    model.factorizations{I}.data = ...
        fmt(model.factorizations{I}.data,true);
    if isstruct(model.factorizations{I}.data)
        cache.factorizations.isincomplete(I) = ...
            model.factorizations{I}.data.incomplete;
        cache.factorizations.issparse(I) = ...
            model.factorizations{I}.data.sparse;
    end
end

% Cache expansion of variables into factors.
cache.factors.expanded = cell(size(model.factors));
cache.factors.sequence = cell(size(model.factors));
cache.factors.state = cell(size(model.factors));
cache.factors.persistent = cell(size(model.factors));
for I = 1:length(model.factors)
    cache.factors.sequence{I} = cell(size(model.factors{I}));
    cache.factors.state{I} = cell(size(model.factors{I}));
    cache.factors.persistent{I} = cell(size(model.factors{I}));
    for J = 1:numel(model.factors{I})
        cache.factors.sequence{I}{J} = cell(1,length(model.factors{I}{J}));
        cache.factors.state{I}{J} = cell(1,length(model.factors{I}{J})-1);
        cache.factors.persistent{I}{J} = cell(1,length(model.factors{I}{J})-1);
    end
end
expand(model.variables,'cache','persistent');

% Cache factors' structure.
cache.factors.structure = ...
    cellfun(@structure,cache.factors.sequence,'UniformOutput',false);

% Cache offsets for variables and expanded factors.
cache.variables.offset = ...
    cumsum([1 cellfun(@(v)numel(serialize(v)),model.variables)]);
cache.factorizations.serialized = cell(size(model.factorizations));
cache.factorizations.offset = cell(size(model.factorizations));
cache.factorizations.suboffset = cell(size(model.factorizations));
for I = 1:length(model.factorizations)
    cache.factorizations.serialized{I} = ...
        serialize(model.factorizations{I}.factors).';
    cache.factorizations.offset{I} = ...
        cumsum([1 cellfun(@(v)numel(serialize(v)), ...
        cache.factors.expanded(cache.factorizations.serialized{I}))]);
    cum = cumsum([0 cellfun(@(f)sum(serialize( ...
        cellfun(@(s)numel(s{end}),f))),cache.factors.sequence( ...
        cache.factorizations.serialized{I}))]);
    fct = cellfun(@(f)cellfun(@(s)false(size(s{end})),f,'UniformOutput',...
        false),cache.factors.sequence( ...
        cache.factorizations.serialized{I}),'UniformOutput',false);
    cache.factorizations.suboffset{I} = cell(1,length(fct));
    for J = 1:length(fct)
        cache.factorizations.suboffset{I}{J} = cell(size(fct{J}));
        for K = 1:numel(fct{J})
            subfct = fct{J};
            subfct{K} = true(size(subfct{K}));
            subfct = find(cell2mat(subfct));
            if all(subfct(2:end) == subfct(1:end-1)+1)
                subfct = [subfct(1) subfct(end)];
            end
            cache.factorizations.suboffset{I}{J}{K} = subfct+cum(J);
        end
    end
end

% Set (relative) weights.
if isempty(options.RelWeights) && isempty(options.Weights)
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
elseif ~isempty(options.RelWeights) && ~isempty(options.Weights)
    warning('sdf:weights',['Both relative and absolute weights ' ...
        'are supplied, proceeding with absolute weights.']);
end
if isempty(options.RelWeights)
    options.RelWeights = ones(1,length(model.factorizations));
end
options.RelWeights = options.RelWeights./sum(options.RelWeights);
if isempty(options.Weights);
    for I = 1:length(model.factorizations)
        if isfield(model.factorizations{I},'data') && ...
           ~strncmpi(model.factorizations{I}.type, 'reg', 3)
            if strcmp(model.factorizations{I}.datatype, 'incomplete')
                NUM = numel(model.factorizations{I}.data.val);
            else 
                NUM = prod(getsize(model.factorizations{I}.data));
            end
        else
            NUM = cache.factorizations.offset{I}(end)-1;
        end
        options.Weights(I) = 2*options.RelWeights(I)/NUM;
    end
end
if length(options.Weights) ~= length(model.factorizations)
    error('sdf:weights',['The number of weights must equal the ' ...
        'number of factorizations.']);
end

% Initialize each factorization's state.
cache.factorizations.state = cell(size(model.factorizations));
for I = 1:length(model.factorizations)
    try 
        model.factorizations{I}.state(I,true);
    catch e
        if strcmpi(e.identifier, 'frob:notImplemented');
            error('sdf_core:notImplemented', ...
                  ['sdf_core does not support the structured tensor type %s, yet. Use ' ...
                   'ful(T) instead in factorization %d.'], type, I);
        end
    end
    if ~strcmpi(options.OptimizationType, 'nls') % if not nls, disable
        model.factorizations{I}.state = @(x,y) 0;
    end
end

% Compute absolute tolerance
if any(strcmp(p.UsingDefaults, 'TolAbs'))
    nrm = 0;
    for I = 1:length(model.factorizations)
        if ~strncmpi(model.factorizations{I}.type, 'reg', 3)
            if model.factorizations{I}.isstructured
                nrm = nrm + 0.5*options.Weights(I)* ...
                      cache.factorizations.state{I}.T2;
            end
        end
    end
    options.TolAbs = nrm * 1e-15;
end

%% Compute bypassing checks
% Compute skipping checks for improved performance
cache.factors.issimple = false(1,length(model.factors));
cache.factors.fastrefs = cell(1,length(model.factors));
for I = 1:length(model.factors)
    cache.factors.issimple(I) = ~any(model.isconstant{I}(:)) && ...
        (length(model.factors{I}) == 1) && ...
        (length(model.factors{I}{1}) == 1);
    if cache.factors.issimple(I)
        seq = model.factors{I}{1}{1};
        dim = cache.factors.structure{I}{1}{1};
        voffset = cache.variables.offset;
        cache.factors.fastrefs{I} = reshape(voffset(seq):voffset(seq+1)-1,dim);
    end
end
% Is simple model? 
cache.bypass.simplemodel = all(cache.factors.issimple);
if cache.bypass.simplemodel
    cache.factors.variablerefs = cellfun(@(f) f{1}{1}, model.factors);
end
% Cache for objective function values
cache.factorizations.objfun = zeros(1, length(model.factorizations));

%% Run optimization algorithm.
if strcmpi(options.OptimizationType, 'nls')
    dF.JHF = @grad;
    dF.JHJx = @JHJx;
    if options.PC && isfield(options,'M') && isa(options.M,'function_handle')
        dF.M = options.M;
    end
else 
    dF = @grad;
    if isfield(options, 'M') && isnan(options.M)
        options = rmfield(options, 'M'); 
    end
end
[z,output] = options.Algorithm(@objfun,dF,model.variables,options);
output.Name = func2str(options.Algorithm);

% Expand the variables into factors.
x = expand(z);

% Compute objective function value.
printwarning = false;
for ii = 1:length(model.factorizations)
    if ~isfield(model.factorizations{ii}, 'data')
        t2 = 0; 
    elseif ~isfield(cache.factorizations.state{ii}, 'T2')
        t2 = frob(model.factorizations{ii}.data);
    else 
        t2 = sqrt(cache.factorizations.state{ii}.T2);
    end
    if strcmp(model.factorizations{ii}.type, 'regL1')
        output.abserr(ii) = 2*cache.factorizations.objfun(ii);        
    elseif strcmp(model.factorizations{ii}.type, 'regL1')
        output.abserr(ii) = cache.factorizations.objfun(ii);        
    else 
        output.abserr(ii) = sqrt(2*cache.factorizations.objfun(ii));
    end
    output.relerr(ii) = output.abserr(ii)/t2;
    
    if isfield(model.factorizations{ii}, 'isstructured') && ...
        model.factorizations{ii}.isstructured && output.relerr(ii) < 5e-8
        printwarning = true;
    end
end
if printwarning
    warning('sdf_core:accuracy', ...
            ['Maximal numerical accuracy for some structured tensors reached. The ' ...
             'result may be improved using ful(T) instead of T as data for ' ...
             'some factorizations.']);
end

% Return output.
if isfield(model.names,'variables')
    sol.variables = cell2struct(z(:),model.names.variables);
else
    sol.variables = z;
end
if isfield(model.names,'factors')
    sol.factors = cell2struct(reshape(expand(z),[],1),model.names.factors);
else
    sol.factors = expand(z);
end
sol.info.Weights = options.Weights;

function x = expand(z,where,evalpersistent)
% Expands the variables into factors stored in the field factors.expanded,
% and stores their computational state in the field factors.state.

    % Save some references for speed.
    factors = model.factors;
    isconst = cache.factors.isconst;
    
    if nargin < 3 || ~ischar(evalpersistent), evalpersistent = false; end
    
    % Expand into cache or into the output variable x.
    if nargin >= 2 && ischar(where)
        
        % For each ith factor...
        for i = 1:length(factors)
            % For each jth subfactor...
            for j = 1:numel(factors{i})
                if isconst{i}(j)
                    % Constant subfactor.
                    cache.factors.sequence{i}{j}{1} = factors{i}{j}{1};
                else
                    % For each kth transformation...
                    cache.factors.sequence{i}{j}{1} = z{factors{i}{j}{1}};
                    for k = 2:length(factors{i}{j})
                        if ~isempty(cache.factors.persistent{i}{j}{k-1})
                            task = struct;
                            task.l = [];
                            task.r = [];
                            task.persistent = ...
                                cache.factors.persistent{i}{j}{k-1};
                        else 
                            task = [];
                        end
                        [cache.factors.sequence{i}{j}{k}, ...
                            cache.factors.state{i}{j}{k-1}] = ...
                            factors{i}{j}{k}( ...
                            cache.factors.sequence{i}{j}{k-1}, task);
                        if evalpersistent
                            if isstruct(cache.factors.state{i}{j}{k-1}) && ...
                                isfield(cache.factors.state{i}{j}{k-1}, ...
                                        'persistent')
                                 cache.factors.persistent{i}{j}{k-1} = ...
                                    cache.factors.state{i}{j}{k-1}.persistent;
                            else
                                cache.factors.persistent{i}{j}{k-1} = [];
                            end
                        end
                    end
                end
            end
            if numel(cache.factors.sequence{i}) == 1
                cache.factors.expanded{i} = ...
                    cache.factors.sequence{i}{j}{end};
            else
                cache.factors.expanded{i} = ...
                    cell2mat(cellfun(@(f)f{end}, ...
                    cache.factors.sequence{i},'UniformOutput',false));
            end
        end
    else
        
        % For each ith factor...
        x = cell(size(factors));
        for i = 1:length(factors)
            % For each jth subfactor...
            sub = cell(size(factors{i}));
            for j = 1:numel(factors{i})
                if isconst{i}(j)
                    % Constant subfactor.
                    sub{j} = factors{i}{j}{1};
                else
                    % For each kth transformation...
                    sub{j} = z{factors{i}{j}{1}};
                    for k = 2:length(factors{i}{j})
                        if ~isempty(cache.factors.persistent{i}{j}{k-1})
                            task = struct;
                            task.l = [];
                            task.r = [];
                            task.persistent = ...
                                cache.factors.persistent{i}{j}{k-1};
                        else 
                            task = [];
                        end                        
                        sub{j} = factors{i}{j}{k}(sub{j},task);
                    end
                end
            end
            if numel(sub) == 1, x{i} = sub{1};
            else x{i} = cell2mat(sub);
            end
        end
        
    end
    
end

function x = derivexpand(r)
% Linearly expand variables using their transformations' Jacobians.

    % Save some references for speed.
    factors = model.factors;
    voffset = cache.variables.offset;
    isconst = cache.factors.isconst;
    sequence = cache.factors.sequence;
    state = cache.factors.state;
    structure = cache.factors.structure;

    % For each ith factor...
    x = cell(size(factors));
    for i = 1:length(factors)
        if cache.factors.issimple(i)
            x{i} = reshape(r(cache.factors.fastrefs{i}), ...
                           size(cache.factors.fastrefs{i}));
            continue;
        end
        % For each jth subfactor...
        sub = cell(size(factors{i}));
        for j = 1:numel(factors{i})
            if isconst{i}(j)
                % Constant subfactor.
                sub{j} = zeros(size(sequence{i}{j}{1}));
            else
                % For each kth transformation...
                seq = factors{i}{j};
                ftr = sequence{i}{j};
                stt = state{i}{j};
                persist = cache.factors.persistent{i}{j};
                dim = structure{i}{j}{1};
                sub{j} = r(voffset(seq{1}):voffset(seq{1}+1)-1);
                if isnumeric(dim)
                    if ~isempty(dim), sub{j} = reshape(sub{j},dim); end
                else sub{j} = deserialize(sub{j},dim);
                end
                for k = 2:length(seq)
                    task = stt{k-1};
                    task.l = [];
                    task.r = sub{j};
                    if ~isempty(persist{k-1})
                        task.persistent = persist{k-1};
                    end
                    sub{j} = seq{k}(ftr{k-1},task);
                end
            end

        end
        if numel(sub) == 1, x{i} = sub{1};
        else x{i} = cell2mat(sub);
        end
    end

end

function y = derivcontract(f,y,l)
% Linearly contract factors using their transformations' Jacobians.

    % Save some references for speed.
    factors = model.factors;
    voffset = cache.variables.offset;
    isconst = cache.factors.isconst;
    sequence = cache.factors.sequence;
    state = cache.factors.state;
    structure = cache.factors.structure;
    serialized = cache.factorizations.serialized{f};
    soffset = cache.factorizations.suboffset{f};

    % Update y with Jacobian-contracted variables.
    for i = 1:length(serialized)
        % For each jth subfactor...
        idx = serialized(i);
        for j = 1:numel(factors{idx})
            
            % Skip this subfactor if it's constant.
            if isconst{idx}(j), continue; end

            % Apply sequence of Jacobian-vector products.
            seq = factors{idx}{j};
            ftr = sequence{idx}{j};
            stt = state{idx}{j};
            persist = cache.factors.persistent{idx}{j};
            if size(soffset{i}{j},2) == 2
                sub = l(soffset{i}{j}(1):soffset{i}{j}(2));
            else
                sub = l(soffset{i}{j});
            end
            if length(seq) > 1
                dim = structure{idx}{j}{end};
                if isnumeric(dim)
                    if ~isempty(dim), sub = reshape(sub,dim); end
                else sub = deserialize(sub,dim);
                end
                for k = length(seq):-1:2
                    task = stt{k-1};
                    task.l = sub;
                    task.r = [];
                    if ~isempty(persist{k-1})
                        task.persistent = persist{k-1};
                    end
                    sub = seq{k}(ftr{k-1},task);
                end
                sub = serialize(sub);
            end
            
            % Update y.
            jdx = voffset(seq{1}):voffset(seq{1}+1)-1;
            y(jdx) = y(jdx)+sub;
        
        end
    end
    
end

function z = getfactors(f,x)
% Retrieves factors of the fth factorization, given the factors x.
    
    % Deserialize the factors.
    if any(strcmp(model.factorizations{f}.type, {'cpd', 'cpdi', 'll1'}))
        fidx = cat(2,model.factorizations{f}.factors{:});
        z = x(fidx);
    else 
        z = deserialize_local(x,model.factorizations{f}.factors);
    end
    function z = deserialize_local(x,dim)
        z = cell(size(dim));
        for i = 1:numel(dim)
            if iscell(dim{i})
                z{i} = deserialize_local(x,dim{i});
            else
                z{i} = x{dim{i}};
            end
        end
    end
    
end

function fval = objfun(z)
    
    % Expand the variables into factors.
    x = expand(z);
    
    % Compute objective function value.
    fval = 0;
    for i = find(options.Weights(:).' ~= 0)
        cache.factorizations.objfun(i) = ...
            model.factorizations{i}.objfun(i,getfactors(i,x));
        fval = fval+options.Weights(i)*cache.factorizations.objfun(i);
    end
end

function grad = grad(z)
    
    % Expand the variables into factors.
    expand(z,'cache');
    
    % Let each factorization save intermediate computations in the cache.
    for i = 1:length(model.factorizations)
        model.factorizations{i}.state(i);
    end
    
    % Compute the gradient.
    grad = zeros(cache.variables.offset(end)-1,1);
    for i = find(options.Weights(:).' ~= 0)
        
        % Compute model's gradient.
        tmp = getfactors(i,cache.factors.expanded);
        JHF = model.factorizations{i}.grad(i,tmp);
        
        % Contract the factor matrices into variables.
        if i == 1
            grad = options.Weights(i)*derivcontract(i,grad,JHF);
        else
            grad = grad+options.Weights(i)*derivcontract(i,grad,JHF);
        end
        
    end
    
end

function y = JHJx(~,x)
    
% Expand the variables into factors.
    if cache.bypass.simplemodel
        x = cellfun(@(idx) reshape(x(idx),size(idx)), cache.factors.fastrefs, ...
                                   'UniformOutput', false);
    else 
        x = derivexpand(x);
    end
    
    % Compute (J(z)'*J(z))*x after expansion.
    y = zeros(cache.variables.offset(end)-1,1);
    for i = find(options.Weights(:).' ~= 0)
        
        % Apply fast matrix-vector product.
        tmpa = getfactors(i,cache.factors.expanded);
        tmpb = getfactors(i,x);
        JHJx = model.factorizations{i}.JHJx(i,tmpa,tmpb);
        
        % Contract the factor matrices into variables.
        if i == 1
            y = options.Weights(i)*derivcontract(i,y,JHJx);
        else
            y = y+options.Weights(i)*derivcontract(i,y,JHJx);
        end
        
    end
    
end

function [z,offset] = deserialize(z,dim,offset)
    if iscell(dim)
        v = z;
        z = cell(size(dim));
        if nargin < 3, offset = 0; end
        for i = 1:numel(z)
            if iscell(dim{i})
                [z{i},offset] = deserialize(v,dim{i},offset);
            else
                n = prod(dim{i}(:));
                z{i} = reshape(v(offset+(1:n)),dim{i});
                offset = offset+n;
            end
        end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
        for i = find(cellfun(@iscell,z(:).'))
            dim{i} = structure(z{i});
        end
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end

% Model: canonical polyadic decomposition ---------------------------------

function state_cpd(f,firstrun)
% Can read from model, cache and options, can save state by writing to the
% structure cache.factorizations.state{f}.

    if nargin == 2 && firstrun
        
        % Store the fraction of known elements.
        if cache.factorizations.isincomplete(f)
            cache.factorizations.scale{f} = ...
                length(model.factorizations{f}.data.val)./...
                prod(model.factorizations{f}.data.size);
        end
        
        % if structure, compute Frobenius norm
        if model.factorizations{f}.isstructured
            cache.factorizations.state{f}.T2 = frob(model.factorizations{f}.data, ...
                                                    'squared');
        end
        
        if strcmpi(model.factorizations{f}.type, 'cpdi')
            [cache.factorizations.j{f}, cache.factorizations.i{f}] = cellfun(@sort, ...
                                                              model.factorizations{f}.data.sub, 'UniformOutput', 0);
            
            cache.factorizations.i{f} = cellfun(@double, cache.factorizations.i{f}, ...
                                                'UniformOutput', false);
            cache.factorizations.j{f} = cellfun(@double, cache.factorizations.j{f}, ...
                                                'UniformOutput', false);
        end
        
        % Set Block-Jacobi preconditioner if
        % - the model is a single CPD and
        % - each factor consists of exactly one subfactor and
        % - no subfactor is transformed and
        % - every subfactor is a reference to a variable and
        % - all variable references are unique and
        % - all factor references are unique and
        % - an NLS algorithm is used.
        options.BJ = false;
        ftr = cell2mat(model.factorizations{f}.factors);
        var = cellfun(@(v)v{1}{1},model.factors(ftr),'UniformOutput',false);
        if options.PC && numel(model.factorizations) == 1 && ...
                all(cellfun(@(v)length(v) == 1,model.factors(ftr))) && ...
                all(cellfun(@(v)length(v{1}) == 1,model.factors(ftr))) && ...
                all(cellfun(@isscalar,var)) && ...
                length(unique(cell2mat(var))) == length(cell2mat(var))&&...
                length(unique(ftr)) == length(ftr) && ...
                ~isempty(strfind(func2str(options.Algorithm),'nls'))
            options.M = @pc_cpd;
            options.BJ = true;
        end
        
        % Skip constant factors
        facts = cat(2, model.factorizations{f}.factors{:});
        constfact = model.isconstant(facts);
        constfact = cellfun(@(c) all(c(:)), constfact);
        cache.factorizations.constfact{f} = constfact;
    end

    % Cache the factor matrices' Gramians.
    idx = cache.factorizations.serialized{f};
    N = length(idx);
    R = size(cache.factors.expanded{idx(1)},2);
    cache.factorizations.state{f}.UHU = zeros(N,R*R);
    for n = 1:N
        tmp = cache.factors.expanded{idx(n)};
        tmp = conj(tmp'*tmp);
        cache.factorizations.state{f}.UHU(n,:) = tmp(:);
    end
    
    % Optionally cache the inverses of the Gramians for the preconditioner.
    % In a faster language, this should be the Cholesky factor instead.
    if options.BJ
        cache.factorizations.state{f}.invW = cell(1,N);
        for n = 1:N
            tmp = cache.factorizations.state{f}.UHU([1:n-1 n+1:N],:);
            if N > 2, tmp = prod(tmp,1); end
            cache.factorizations.state{f}.invW{n} = ...
                inv(reshape(tmp,[R R]));
        end
    end
    
    if strcmpi(model.factorizations{f}.type, 'cpdi')
        J = cell(R,N);
        T = model.factorizations{f}.data;
        U = cache.factors.expanded(idx);
        for n = 1:N
            dims = [1:n-1 n+1:N];
            i = cache.factorizations.i{f}{n};
            j = cache.factorizations.j{f}{n};
            for r = 1:R
                val = U{dims(1)}(T.sub{dims(1)},r);
                for k = dims(2:end)
                    val = val .* U{k}(T.sub{k},r);
                end
                val = val(i);
                tmp = sparse(i,j,val,length(T.val),T.size(n));
                J{r,n} = tmp;                
            end
        end
        cache.factorizations.state{f}.J = cat(2, J{:});
    end 
    
end

function fval = objfun_cpd(f,z)
    
    % CPD objective function.
    T = model.factorizations{f}.data;
    isincomplete = cache.factorizations.isincomplete(f);
    issparse = cache.factorizations.issparse(f);
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = z{1}*kr(z(end:-1:2)).';
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:size(z{1},2)
            tmp = z{1}(T.sub{1},r);
            for n = 2:length(z), tmp = tmp.*z{n}(T.sub{n},r); end
            fval = fval+tmp;
        end
    end
    if isincomplete
        T.val = fval;
        if ~isempty(T.matrix)
            T.matrix = sparse(double(T.sub{1}), ...
                double(1+idivide(T.ind-1,int64(size(T.matrix,1)))), ...
                double(fval),size(T.matrix,1),size(T.matrix,2));
        end
        cache.factorizations.residual{f} = T;
    else
        if issparse, size_tens = T.size;
        else size_tens = size(T); end
        cache.factorizations.residual{f} = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function fval = objfun_cpd_structured(f,z)
    T = model.factorizations{f}.data;
    if isfield(cache.factorizations.state{f}, 'T2')
        T2 = cache.factorizations.state{f}.T2;
    else 
        T2 = frob(T, 'squared');
        cache.factorizations.state{f}.T2 = T2;
    end
    fval = abs(0.5*T2 - real(inprod(T,z)) + 0.5*frob(z, 'squared'));
end

function grad = grad_cpd(f,z)
    
    % CPD scaled conjugate cogradient.
    E = cache.factorizations.residual{f};
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    grad = zeros(offset(end)-1,1);
    for n = 1:length(z)
        if ~constfact(n)
            tmp = full(mtkrprod(E,z,n));
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end
end

function grad = grad_cpd_structured(f,z)
    
    % CPD scaled conjugate cogradient.
    T = model.factorizations{f}.data;
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    grad = zeros(offset(end)-1,1);
    for n = 1:length(z)
        if ~constfact(n)
            tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end
    
end

function y = JHJx_cpd(f,z,x)
    
    % CPD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    R = size(z{1},2);
    N = length(z);
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    UHU = cache.factorizations.state{f}.UHU;
    XHU = zeros(R,R,N);
    y = zeros(offset(end)-1,1);
    for n = 1:N
        if ~constfact(n)
            Wn = UHU([1:n-1 n+1:N],:);
            if N > 2, Wn = prod(Wn,1); end
            XHU(:,:,n) = conj(x{n}'*z{n});
            y(offset(n):offset(n+1)-1) = x{n}*reshape(Wn,[R R]);
        end
    end
    for n = 1:N-1
        if ~constfact(n)
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
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(f)
        y = y*cache.factorizations.scale{f};
    end
end

function y = JHJxi_cpd(f,~,x)
    J = cache.factorizations.state{f}.J;
    x = cellfun(@(x) x(:), x, 'UniformOutput', false);
    x = cat(1, x{:});
    y = J*x;
    y = conj(y'*J);
    y = y(:);
end

function x = pc_cpd(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = zeros(size(b));
    offset = cache.factorizations.offset{1};
    invW = cache.factorizations.state{1}.invW;
    for n = 1:length(offset)-1
        idx = offset(n):offset(n+1)-1;
        tmp = reshape(b(idx),[],size(invW{1},1))*invW{n};
        x(idx) = tmp(:);
    end
    x = x/options.Weights;
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(1)
        x = x/cache.factorizations.scale{1};
    end
    
end

% Model: canonical polyadic decomposition ---------------------------------

function state_ll1(f,firstrun)
% Can read from model, cache and options, can save state by writing to the
% structure cache.factorizations.state{f}.

    idx = cache.factorizations.serialized{f};
    N = length(idx);
    
    if nargin == 2 && firstrun
        
        % Store the fraction of known elements.
        if cache.factorizations.isincomplete(f)
            cache.factorizations.scale{f} = ...
                length(model.factorizations{f}.data.val)./...
                prod(model.factorizations{f}.data.size);
        end
        
        % if structure, compute Frobenius norm
        if model.factorizations{f}.isstructured
            cache.factorizations.state{f}.T2 = frob(model.factorizations{f}.data, ...
                                                    'squared');
        end
        
        % Set Block-Jacobi preconditioner if
        % - the model is a single LL1 and
        % - each factor consists of exactly one subfactor and
        % - no subfactor is transformed and
        % - every subfactor is a reference to a variable and
        % - all variable references are unique and
        % - all factor references are unique and
        % - an NLS algorithm is used.
        options.BJ = false;
        ftr = cell2mat(model.factorizations{f}.factors);
        var = cellfun(@(v)v{1}{1},model.factors(ftr),'UniformOutput',false);
        if options.PC && numel(model.factorizations) == 1 && ...
                all(cellfun(@(v)length(v) == 1,model.factors(ftr))) && ...
                all(cellfun(@(v)length(v{1}) == 1,model.factors(ftr))) && ...
                all(cellfun(@isscalar,var)) && ...
                length(unique(cell2mat(var))) == length(cell2mat(var))&&...
                length(unique(ftr)) == length(ftr) && ...
                ~isempty(strfind(func2str(options.Algorithm),'nls'))
            options.M = @pc_ll1;
            options.BJ = true;
        end
        
        % Compute expansion vectors and contraction P
        L = model.factorizations{f}.L;
        expvec = arrayfun(@(n) ones(1, L(n))*n, 1:length(L), 'UniformOutput', false);
        expvec = cat(2, expvec{:});
        cache.factorizations.state{f}.expvec = expvec;
        cache.factorizations.state{f}.P = sparse(1:sum(L), expvec, 1);
        
        % Skip constant factors
        facts = cat(2, model.factorizations{f}.factors{:});
        constfact = model.isconstant(facts);
        constfact = cellfun(@(c) all(c(:)), constfact);
        cache.factorizations.constfact{f} = constfact;
    end

    % Cache the factor matrices' Gramians.
    expvec = cache.factorizations.state{f}.expvec;
    P      = cache.factorizations.state{f}.P;
    z      = cache.factors.expanded(idx);
    z{N} = z{N}(:,expvec);
    UHU = cellfun(@(u) conj(u'*u), z, 'UniformOutput', false);
    cache.factorizations.state{f}.UHU = UHU;
    
    % Optionally cache the inverses of the Gramians for the preconditioner.
    % In a faster language, this should be the Cholesky factor instead.
    if options.BJ
        cache.factorizations.state{f}.invW = cell(1,N);
        cache.factorizations.state{f}.invW{1} = inv(UHU{2}.*UHU{3});
        cache.factorizations.state{f}.invW{2} = inv(UHU{1}.*UHU{3});
        cache.factorizations.state{f}.invW{3} = P'*inv(UHU{1}.*UHU{2}).'*P;
    end
end

function fval = objfun_ll1(f,z)
% LL1 objective function.
    N = length(z);
    T = model.factorizations{f}.data;
    z{N} = z{N}(:, cache.factorizations.state{f}.expvec);
    cache.factorizations.residual{f} = cpdres(T,z,struct('Format', false));
    if isstruct(cache.factorizations.residual{f}), 
        fval = cache.factorizations.residual{f}.val;
    else 
        fval = cache.factorizations.residual{f};
    end
    fval = 0.5*(fval(:)'*fval(:));
end

function fval = objfun_ll1_structured(f,z)
    N = length(z);
    T = model.factorizations{f}.data;
    z{N} = z{N}(:, cache.factorizations.state{f}.expvec);
    if isfield(cache.factorizations.state{f}, 'T2')
        T2 = cache.factorizations.state{f}.T2;
    else 
        T2 = frob(T, 'squared');
        cache.factorizations.state{f}.T2 = T2;
    end
    fval = abs(0.5*T2 - real(inprod(T,z)) + 0.5*frob(z, 'squared'));
end

function grad = grad_ll1(f,z)
% LL1 scaled conjugate cogradient.
    N = length(z);
    E = cache.factorizations.residual{f};
    offset = cache.factorizations.offset{f};
    grad = nan(offset(end)-1, 1);
    z{N} = z{N}(:, cache.factorizations.state{f}.expvec);
    P = cache.factorizations.state{f}.P;
    constfact = cache.factorizations.constfact{f};
    for n = 1:length(z)
        if ~constfact(n)
            tmp = full(mtkrprod(E, z, n));
            if n == length(z), tmp = tmp*P; end
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end  
end

function grad = grad_ll1_structured(f,z)
% LL1 scaled conjugate cogradient.
    N = length(z);
    offset = cache.factorizations.offset{f};
    T = model.factorizations{f}.data;
    grad = nan(offset(end)-1, 1);
    z{N} = z{N}(:, cache.factorizations.state{f}.expvec);
    P = cache.factorizations.state{f}.P;
    constfact = cache.factorizations.constfact{f};
    for n = 1:length(z)
        if ~constfact(n)
            tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            if n == length(z), tmp = tmp*P; end
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end  
end

function y = JHJx_ll1(f,z,x)
% LL1 fast Jacobian's Gramian vector product.
% Ignores the fact that the tensor might be incomplete.
    N = length(z);
    constfact = cache.factorizations.constfact{f};
    offset = cache.factorizations.offset{f};
    expvec = cache.factorizations.state{f}.expvec;
    P      = cache.factorizations.state{f}.P;
    UHU    = cache.factorizations.state{f}.UHU;
    XHU    = cell(1,N);
    y      = zeros(offset(end)-1,1);
    L      = model.factorizations{f}.L;

    % Diagonal blocks
    for n = 1:3
        if ~constfact(n)
            idxw = true(1,N);
            idxw(n) = false;
            Wn = prod(cat(3,UHU{idxw}),3);
            if n == 3, Wn = P'*Wn*P; end
            tmp = x{n};
            XHU{n} = conj(tmp'*z{n});
            y(offset(n):offset(n+1)-1) = tmp*Wn;
        end
    end
    XHU{3} = XHU{3}(expvec,expvec);
    
    % Off diagonal blocks
    for n = 1:N-1
        if ~constfact(n)
            idxn = offset(n):offset(n+1)-1;
            Wn = zeros(sum(L));
            for m = n+1:N
                idxm = offset(m):offset(m+1)-1;
                idxw = true(1,N);
                idxw([n m]) = false;
                Wn = Wn + UHU{idxw}.*XHU{m};
                Wnm = UHU{idxw}.*XHU{n};
                if m == 3, Wnm = P'*Wnm*P; end
                tmp = z{m}*Wnm;
                y(idxm) = y(idxm) + tmp(:);
            end
            tmp = z{n}*Wn;
            y(idxn) = y(idxn) + tmp(:);
        end
    end

    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(f)
        y = y*cache.factorizations.scale{f};
    end
    
end

function x = pc_ll1(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = zeros(size(b));
    offset = cache.factorizations.offset{1};
    invW = cache.factorizations.state{1}.invW;
    for n = 1:length(offset)-1
        idx = offset(n):offset(n+1)-1;
        tmp = reshape(b(idx),[],size(invW{n},1))*invW{n};
        x(idx) = tmp(:);
    end
    x = x/options.Weights;
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(1)
        x = x/cache.factorizations.scale{1};
    end
    
end

% Model: block term decomposition -----------------------------------------

function state_btd(f,firstrun)
% Can read from model, cache and options, can save state by writing to the
% structure cache.factorizations.state{f}.
    
    if nargin == 2 && firstrun
        
        % Store the fraction of known elements.
        if cache.factorizations.isincomplete(f)
            cache.factorizations.scale{f} = ...
                length(model.factorizations{f}.data.val)./...
                prod(model.factorizations{f}.data.size);
        end
        
        % if structure, compute Frobenius norm
        if model.factorizations{f}.isstructured
            cache.factorizations.state{f}.T2 = frob(model.factorizations{f}.data, ...
                                                    'squared');
        end
        
        % Set Block-Jacobi preconditioner if
        % - the model is a single BTD and
        % - each factor consists of exactly one subfactor and
        % - no subfactor is transformed and
        % - every subfactor is a reference to a variable and
        % - all variable references are unique and
        % - all factor references are unique and
        % - an NLS algorithm is used.
        options.BJ = false;
        ftr = cache.factorizations.serialized{1};
        var = cellfun(@(v)v{1}{1},model.factors(ftr),'UniformOutput',false);
        if options.PC && numel(model.factorizations) == 1 && ...
                all(cellfun(@(v)length(v) == 1,model.factors(ftr))) && ...
                all(cellfun(@(v)length(v{1}) == 1,model.factors(ftr))) && ...
                all(cellfun(@isscalar,var)) && ...
                length(unique(cell2mat(var))) == length(cell2mat(var))&&...
                length(unique(ftr)) == length(ftr) && ...
                ~isempty(strfind(func2str(options.Algorithm),'nls'))
            options.M = @pc_btd;
            options.BJ = true;
        end
        
        % Skip constant factors
        facts = cat(2, model.factorizations{f}.factors{:});
        facts = cat(2, facts{:});
        constfact = model.isconstant(facts);
        constfact = cellfun(@(c) all(c(:)), constfact);
        cache.factorizations.constfact{f} = constfact;
    end
    
    % Cache the factor matrices' Gramians.
    U = getfactors(f,cache.factors.expanded);
    R = length(U);
    N = length(U{1})-1;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    cache.factorizations.state{f}.UHU = ...
        arrayfun(@(i,n,j)U{i}{n}'*U{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    
    % Optionally cache some results for the block-Jacobi preconditioner.
    if options.BJ
        [idx,jdx] = ndgrid(1:R,1:N);
        UHU = cache.factorizations.state{f}.UHU;
        cache.factorizations.state{f}.invSKS = arrayfun( ...
            @(r,n)inv(mtkronprod(U{r}{end},UHU(r,:,r),n)* ...
            conj(reshape(permute(U{r}{end},[1:n-1 n+1:N n]), ...
            [],size(U{r}{end},n)))),idx,jdx,'UniformOutput',false);
        cache.factorizations.state{f}.invUHU = arrayfun( ...
            @(r,n)inv(UHU{r,n,r}),idx,jdx,'UniformOutput',false);
    end

end

function fval = objfun_btd(f,z)
    
    % BTD objective function.
    T = model.factorizations{f}.data;
    isincomplete = cache.factorizations.isincomplete(f);
    issparse = cache.factorizations.issparse(f);
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = z{1}{1}*mtkronprod(z{1}{end},z{1}(1:end-1),1,'H');
        for r = 2:length(z)
            fval = fval+z{r}{1}*mtkronprod(z{r}{end},z{r}(1:end-1),1,'H');
        end
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:length(z)
            size_core = cellfun('size',z{r}(1:end-1),2);
            idx = cell(1,length(size_core));
            S = z{r}{end};
            for i = 1:numel(S)
                [idx{:}] = ind2sub(size_core,i);
                tmp = S(idx{:})*z{r}{1}(T.sub{1},idx{1});
                for n = 2:length(size_core)
                    tmp = tmp.*z{r}{n}(T.sub{n},idx{n});
                end
                fval = fval+tmp;
            end
        end
    end
    if isincomplete
        T.val = fval;
        if ~isempty(T.matrix)
            T.matrix = sparse(double(T.sub{1}), ...
                double(1+idivide(T.ind-1,int64(size(T.matrix,1)))), ...
                double(fval),size(T.matrix,1),size(T.matrix,2));
        end
        cache.factorizations.residual{f} = T;
    else
        if issparse, size_tens = T.size;
        else size_tens = size(T); end
        cache.factorizations.residual{f} = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function fval = objfun_btd_structured(f,z)
    T = model.factorizations{f}.data;
    if isfield(cache.factorizations.state{f}, 'T2')
        T2 = cache.factorizations.state{f}.T2;
    else 
        T2 = frob(T, 'squared');
        cache.factorizations.state{f}.T2 = T2;
    end
    fval = abs(0.5*T2 - real(inprod(T,z)) + 0.5*frob(z, 'squared'));
end

function grad = grad_btd(f,z)
    
    % BTD scaled conjugate cogradient.
    N = length(z{1})-1;
    E = cache.factorizations.residual{f};
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    grad = zeros(offset(end)-1,1);
    cnt = 1;
    for r = 1:length(z)
        U = z{r}(1:N);
        S = conj(z{r}{end});
        for n = 1:N
            if ~constfact(cnt)
                tmp = full(mtkronprod(E,U,n))* ...
                      reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));
                grad(offset(cnt):offset(cnt+1)-1) = tmp(:);
            end
            cnt = cnt+1;
        end
        if ~constfact(cnt)
            tmp = full(mtkronprod(E,U,0));
            grad(offset(cnt):offset(cnt+1)-1) = tmp;
        end
        cnt = cnt+1;
    end
    
end

function grad = grad_btd_structured(f,z)
% BTD scaled conjugate cogradient.
    N = length(z{1})-1;
    T = model.factorizations{f}.data;
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    grad = zeros(offset(end)-1,1);
    cnt = 1;
    for r = 1:length(z)
        U = z{r}(1:N);
        S = conj(z{r}{end});
        for n = 1:N
            if ~constfact(cnt)
                tmp = full(mtkronprod(z,U,n)-mtkronprod(T,U,n))* ...
                      reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));
                grad(offset(cnt):offset(cnt+1)-1) = tmp(:);
            end
            cnt = cnt+1;
        end
        if ~constfact(cnt)
            tmp = full(mtkronprod(z,U,0)-mtkronprod(T,U,0));
            grad(offset(cnt):offset(cnt+1)-1) = tmp;
        end
        cnt = cnt+1;
    end
end

function y = JHJx_btd(f,z,x)
    
    % BTD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    R = length(z);
    N = length(z{1})-1;
    offset = cache.factorizations.offset{f};
    constfact = cache.factorizations.constfact{f};
    UHU = cache.factorizations.state{f}.UHU;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    XHU = arrayfun(@(i,n,j)x{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    y = zeros(offset(end)-1,1);
    cnt = 1;
    for ri = 1:R
        % Factor matrices.
        for ni = 1:N
            if ~constfact(cnt) 
                idx = offset(cnt):offset(cnt+1)-1;
                Sri = permute(z{ri}{end},[1:ni-1 ni+1:N ni]);
                Sri = conj(reshape(Sri,[],size(Sri,N)));
                for rj = 1:R
                    Srj = z{rj}{end};
                    tmp = mtkronprod(x{rj}{end},UHU(rj,:,ri),ni);
                    for nj = [1:ni-1 ni+1:N]
                        proj = UHU(rj,:,ri);
                        proj{nj} = XHU{rj,nj,ri};
                        tmp = tmp+mtkronprod(Srj,proj,ni);
                    end
                    tmp = z{rj}{ni}*(tmp*Sri);
                    tmp = tmp+x{rj}{ni}* ...
                          (mtkronprod(z{rj}{end},UHU(rj,:,ri),ni)*Sri);
                    if rj == 1, y(idx) = tmp(:);
                    else y(idx) = y(idx)+tmp(:); end
                end
            end
            cnt = cnt+1;
        end
        
        % Core tensor.
        if ~constfact(cnt)
            idx = offset(cnt):offset(cnt+1)-1;
            for rj = 1:R
                Srj = z{rj}{end};
                tmp = mtkronprod(x{rj}{end},UHU(rj,:,ri),0);
                for nj = 1:N
                    proj = UHU(rj,:,ri);
                    proj{nj} = XHU{rj,nj,ri};
                    tmp = tmp+mtkronprod(Srj,proj,0);
                end
                if rj == 1, y(idx) = tmp(:);
                else y(idx) = y(idx)+tmp(:); end
            end
        end
        cnt = cnt+1;
    end

    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(f)
        y = y*cache.factorizations.scale{f};
    end
    
end

function x = pc_btd(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    x = zeros(size(b));
    R = length(model.factorizations{1}.factors);
    N = length(model.factorizations{1}.factors{1})-1;
    offset = cache.factorizations.offset{1};
    invSKS = cache.factorizations.state{1}.invSKS;
    invUHU = cache.factorizations.state{1}.invUHU;
    cnt = 1;
    for r = 1:R
        for n = 1:N
            idx = offset(cnt):offset(cnt+1)-1;
            tmp = reshape(b(idx),[],size(invSKS{r,n},1))*invSKS{r,n};
            x(idx) = tmp(:);
            cnt = cnt+1;
        end
        idx = offset(cnt):offset(cnt+1)-1;
        size_core = cellfun('size',invUHU(r,:),1);
        x(idx) = mtkronprod(reshape(b(idx),size_core),invUHU(r,:),0);
        cnt = cnt+1;
    end
    x = x/options.Weights;
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(1)
        x = x/cache.factorizations.scale{1};
    end
    
end

% Model: L2 regularization ------------------------------------------------

function state_regL2(f,firstrun)
    
    % Format right hand side, if available.
    if nargin == 2 && firstrun
         cache.factorizations.hasdata(f) = ...
             isfield(model.factorizations{f},'data');
         if cache.factorizations.hasdata(f) && ...
            ~iscell(model.factorizations{f}.data)
             model.factorizations{f}.data = ...
                 {model.factorizations{f}.data};
         end
         if cache.factorizations.hasdata(f)
             model.factorizations{f}.data = cellfun(@(d)full(d), ...
                 model.factorizations{f}.data,'UniformOutput',false);
         end
    end
    
end

function fval = objfun_regL2(f,z)
    % L2 regularization objective function.
    fval = 0;
    for i = 1:length(z)
        e = z{i}(:);
        if cache.factorizations.hasdata(f)
            e = e-model.factorizations{f}.data{i}(:);
        end
        fval = fval+(e'*e);
    end
    fval = 0.5*fval;
end

function grad = grad_regL2(f,z)
    % L2 regularization scaled conjugate cogradient.
    if cache.factorizations.hasdata(f)
        grad = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        grad = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
end

function y = JHJx_regL2(~,~,x)
    % L2 regularization Hessian vector product.
    y = cell2mat(cellfun(@(f)f(:),x(:),'UniformOutput',false));
end

% Model: L1 regularization ------------------------------------------------

function state_regL1(f,firstrun)
    
    % Format right hand side, if available.
    if nargin == 2 && firstrun
         cache.factorizations.hasdata(f) = ...
             isfield(model.factorizations{f},'data');
         if cache.factorizations.hasdata(f) && ...
            ~iscell(model.factorizations{f}.data)
             model.factorizations{f}.data = ...
                 {model.factorizations{f}.data};
         end
         if cache.factorizations.hasdata(f)
             model.factorizations{f}.data = cellfun(@(d)full(d), ...
                 model.factorizations{f}.data,'UniformOutput',false);
         end
         if ~isfield(options,'mu')
             if cache.factorizations.hasdata(f)
                m = max(abs(serialize(model.factorizations{f}.data)));
             else
                m = 0;
             end
             options.Mu = max(m,1)/100;
         end
    end
    
end

function fval = objfun_regL1(f,z)
    % Approximate L1 regularization objective function.
    fval = 0;
    for i = 1:length(z)
        e = z{i}(:);
        if cache.factorizations.hasdata(f)
            e = e-model.factorizations{f}.data{i}(:);
        end
        far = abs(e) > options.Mu;
        e(far) = abs(e(far));
        e2 = e(~far);
        e2 = e2.*conj(e2);
        e(~far) = 2*options.Mu*e2./(e2+options.Mu^2);
        fval = fval+sum(e);
    end
    fval = 0.5*fval;
end

function grad = grad_regL1(f,z)
    % Approximate L1 regularization scaled conjugate cogradient.
    if cache.factorizations.hasdata(f)
        grad = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        grad = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    far = abs(grad) > options.Mu;
    grad(far) = grad(far)./(2*abs(grad(far)));
    grad(~far) = 2*options.Mu^3*grad(~far)./ ...
        (grad(~far).*conj(grad(~far))+options.Mu^2).^2;
end

function y = JHJx_regL1(f,z,x)
    % Approximate L1 regularization Hessian vector product.
    if cache.factorizations.hasdata(f)
        hess = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        hess = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    if ~isreal(hess)
        error('sdf_nls:regL1',['Please use sdf_minf when applying L1 ' ...
            'regularization on complex factors.']);
    end
    far = abs(hess) > options.Mu;
    hess(far) = 0;
    if any(hess)
        hess(~far) = 2*options.Mu^3*(options.Mu^2-3*hess(~far).^2)./ ...
            (options.Mu^2+hess(~far).^2).^3;
        y = hess.*cell2mat(cellfun(@(f)f(:),x(:),'UniformOutput',false));
    else
        y = hess;
    end
end


% Model: L0 regularization ------------------------------------------------

function state_regL0(f,firstrun)
    % Format right hand side, if available.
    if nargin == 2 && firstrun
         cache.factorizations.hasdata(f) = ...
             isfield(model.factorizations{f},'data');
         if cache.factorizations.hasdata(f) && ...
            ~iscell(model.factorizations{f}.data)
             model.factorizations{f}.data = ...
                 {model.factorizations{f}.data};
         end
         if cache.factorizations.hasdata(f)
             model.factorizations{f}.data = cellfun(@(d)full(d), ...
                 model.factorizations{f}.data,'UniformOutput',false);
         end
         if isfield(model.factorizations{f}, 'sigma')
             cache.factorizations.state{f}.sigma = model.factorizations{f}.sigma;
         else 
             cache.factorizations.state{f}.sigma = 0.01;
         end
         if isnumeric(cache.factorizations.state{f}.sigma)
             cache.factorizations.state{f}.sigma = @(k,z) cache.factorizations.state{f}.sigma;
         end
         cache.factorizations.state{f}.iteration = -1;
    end
    cache.factorizations.state{f}.iteration = ...
        cache.factorizations.state{f}.iteration + 1;
end

function fval = objfun_regL0(f,z)
    % Approximate L0 regularization objective function.
    fval = 0;
    k = cache.factorizations.state{f}.iteration;
    sigma = cache.factorizations.state{f}.sigma(k,z);
    for i = 1:length(z)
        e = z{i}(:);
        if cache.factorizations.hasdata(f)
            e = e-model.factorizations{f}.data{i}(:);
        end
        e = 1 - exp(-e.^2/sigma^2);
        fval = fval + sum(e(:));
    end
end

function grad = grad_regL0(f,z)
    % Approximate L0 regularization scaled conjugate cogradient.
    k = cache.factorizations.state{f}.iteration;
    sigma = cache.factorizations.state{f}.sigma(k,z);
    if cache.factorizations.hasdata(f)
        grad = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        grad = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    grad = exp(-grad.^2/sigma^2)*2.*grad/sigma^2;
end

function y = JHJx_regL0(f,z,x)
% Approximate L0 regularization Hessian vector product.
    k = cache.factorizations.state{f}.iteration;
    sigma = cache.factorizations.state{f}.sigma(k,z);
    if cache.factorizations.hasdata(f)
        hess = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        hess = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    if ~isreal(hess)
        error('sdf_nls:regL0',['Please use sdf_minf when applying L0 ' ...
            'regularization on complex factors.']);
    end
    hess = -2/sigma^2*exp(-hess.^2/sigma^2).*(1 - 2/sigma^2*hess.^2);
    y = hess.*cell2mat(cellfun(@(f)f(:),x(:),'UniformOutput',false));
end

end
