function [varargout] = sdf_check(model, varargin)
%SDF_CHECK SDF Language parser and syntax/consistency checker.
%   ISCORRECT=SDF_CHECK(MODEL) returns true if the given model is
%   syntactically correct, and if the given model is consistent. A model is
%   consistent if all variables and factors references can be resolved, all
%   factors can be constructed (no transformation errors, and no
%   concatenation errors), and all factors in a factorization are compatible
%   with the factorization and, if provided, the data. If the model is
%   correct (syntax and consistency), the optimization routines can handle
%   the model. 
%
%   SDF_CHECK(MODEL, 'print') prints the interpretated model, and allows the
%   user to see all conversions made from variable to factors, and from
%   factors to factorizations.
%   
%   The remaining options are only relevant for algorithm developers.
%
%   OUT = SDF_CHECK(MODEL, 'internal') translates the given model in the
%   external DSL (available to the user) to the more restricted internal DSL 
%   (available for algorithm developers). If no errors are thrown, OUT is a
%   valid model.
%
%   OUT = SDF_CHECK(MODEL, 'internal', restriction1, restriction2, ...) can
%   be used to disable parts of the syntax. The following string restrictions
%   are defined: 
%      - 'noTransformations'       no transformations allowed
%      - 'noCPD'                   no CPD or CPDI allowed as factorization
%      - 'noLL1'                   no LL1 allowed as factorization
%      - 'noBTD'                   no BTD allowed as factorization
%      - 'noLMLRA'                 no LMLRA allowed as factorization
%      - 'noRegL0'                 no regL0 allowed as factorization
%      - 'noRegL1'                 no regL1 allowed as factorization
%      - 'noRegL2'                 no regL2 allowed as factorization
%      - 'noRegularization'        no regL0, regL1 or regL2 allowed as factorization
%      - 'onlyCPD'                 only CPD and CPDI allowed as factorization
%      - 'onlyLL1'                 only LL1 allowed as factorization
%      - 'onlyBTD'                 only BTD allowed as factorization
%      - 'onlyLMLRA'               only LMLRA allowed as factorization
%      - 'onlySimpleFactors'       no transformations or concatenations allowed

%   Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%              Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version History:
%   - 2015/12/01   NV      Initial version

%% Process options
% add values for string options without values
stringoptions = {'print', 'internal', 'noTransformations', 'noCPD', 'noCPDI', ...
                 'noBTD', 'noLL1', 'noRegularization', 'noRegL0', 'noRegL1', ...
                 'noRegL2', 'onlyCPD', 'onlyLL1', 'onlyBTD', 'onlySimpleFactors', ...
                 'expand'};
k = 1;
while k <= length(varargin)
    idx = cellfun(@(s) strcmpi(varargin{k}, s), stringoptions);
    if any(idx) && (k == length(varargin) || ~islogical(varargin{k+1}));
        varargin = [varargin(1:k), true, varargin(k+1:end)];
        k = k + 2;
    else 
        k = k + 1;
    end 
end 

% Parse options
p = inputParser;
p.addOptional('StopOnError', true);
p.addOptional('Print', false);
p.addOptional('PrintErrors', true);
p.addOptional('Internal', false);
p.addOptional('Expand', false);
p.addOptional('NoTransformations', false); % disallow transformations
p.addOptional('OnlySimpleFactors', false); % no transf, no concat.
p.addOptional('NoCPD', false); 
p.addOptional('NoCPDI', false); 
p.addOptional('NoLL1', false); 
p.addOptional('NoBTD', false); 
p.addOptional('NoLMLRA', false); 
p.addOptional('NoRegularization', false); 
p.addOptional('NoRegL0', false); 
p.addOptional('NoRegL1', false); 
p.addOptional('NoRegL2', false); 
p.addOptional('OnlyCPD', false); 
p.addOptional('OnlyLL1', false); 
p.addOptional('OnlyBTD', false); 
p.addOptional('OnlyLMLRA', false); 
p.parse(varargin{:});
options = p.Results;

if options.Internal && options.Expand
    error('sdf_check:incompatibleOptions', ['The options Internal and Expand ' ...
                        'cannot be true at the same time']);
end

errorlist = {};
isfunc = @(f) isa(f, 'function_handle');
findname = @(name, list) find(cellfun(@(l) strcmp(name, l), list));
knownfactorizations = {'cpd', 'cpdi', 'll1', 'btd', 'lmlra', 'regL0', 'regL1', 'regL2'};
allowedfactorizations = {};
if options.OnlyCPD && options.NoCPD
    error('sdf_check:incompatibleOptions', ['The options OnlyCPD and noCPD ' ...
                        'cannot be combined']);
elseif options.OnlyCPD
    options.NoBTD = true;
    options.NoLMLRA = true;
    options.NoLL1 = true;
    options.NoRegL0 = true;
    options.NoRegL1 = true;
    options.NoRegL2 = true;
end
if options.OnlyLL1 && options.NoLL1
    error('sdf_check:incompatibleOptions', ['The options OnlyLL1 and noLL1 ' ...
                        'cannot be combined']);
elseif options.OnlyLL1
    options.NoCPD = true;
    options.NoBTD = true;
    options.NoLMLRA = true;
    options.NoRegL0 = true;
    options.NoRegL1 = true;
    options.NoRegL2 = true;
end
if options.OnlyBTD && options.NoBTD
    error('sdf_check:incompatibleOptions', ['The options OnlyBTD and noBTD ' ...
                        'cannot be combined']);
elseif options.OnlyBTD
    options.NoCPD = true;
    options.NoLMLRA = true;
    options.NoLL1 = true;
    options.NoRegL0 = true;
    options.NoRegL1 = true;
    options.NoRegL2 = true;
end
if options.OnlyLMLRA && options.NoLMLRA
    error('sdf_check:incompatibleOptions', ['The options OnlyLMLRA and noLMLRA ' ...
                        'cannot be combined']);
elseif options.OnlyLMLRA
    options.NoCPD = true;
    options.NoBTD = true;
    options.NoLL1 = true;
    options.NoRegL0 = true;
    options.NoRegL1 = true;
    options.NoRegL2 = true;
end
if ~options.NoCPD, 
    allowedfactorizations{end+1} = 'cpd'; 
    if ~options.NoCPDI
        allowedfactorizations{end+1} = 'cpdi'; 
    end
end
if ~options.NoLL1, allowedfactorizations{end+1} = 'll1'; end
if ~options.NoBTD,
    allowedfactorizations{end+1} = 'btd'; 
end
if ~options.NoLMLRA,
    allowedfactorizations{end+1} = 'lmlra'; 
end
if ~options.NoRegL0 && ~options.NoRegularization, 
    allowedfactorizations{end+1} = 'regL0'; 
end
if ~options.NoRegL1 && ~options.NoRegularization, 
    allowedfactorizations{end+1} = 'regL1'; 
end
if ~options.NoRegL2 && ~options.NoRegularization, 
    allowedfactorizations{end+1} = 'regL2'; 
end
if options.OnlySimpleFactors, options.NoTransformations = true; end

reservedKeywords = {'(constant)'};
modelname = inputname(1);
isconsistent = true;

[~,~,~,~,v] = regexp(version('-release'),'([0-9]+)([ab])');
displaylinks = (usejava('Desktop') && str2double(v{1}{1}) > 2011) || ...
    (str2double(v{1}{1}) == 2011 && strcmpi(v{1}{2},'b'));

%% Start model check
try   
    %% Valid model
    if ~isstruct(model),
        adderror('sdf:syntax:invalidModel', 'The model should be a struct');
    end
    
    %% Options
    if exist('tensorlabsettings.mat', 'file')
        % if a settings file can be found, apply sdf options if set
        settings = load('tensorlabsettings.mat');
        if isfield(settings, 'sdfoptions')
            if ~isfield(model,'options') || ~isstruct(model.options)
                model.options = settings.sdfoptions;
            else
                % merge given options and settings, overwrite settings if given
                fn1 = fieldnames(model.options);
                fn2 = fieldnames(settings.sdfoptions);
                mask = cellfun(@(s) ~any(strcmpi(fn1, s)), fn2);
                data1 = struct2cell(model.options);
                data2 = struct2cell(settings.sdfoptions);
                model.options = cell2struct([data1; data2(mask)], [fn1; fn2(mask)]);
            end
        end 
    end
    if ~isfield(model, 'options')
        model.options = struct;
    end
    if ~isfield(model.options, 'ImplicitFactors')
        model.options.ImplicitFactors = false;
    end
    
    %% Completeness
    fields = fieldnames(model);
    if ~ismember(fields, 'variables')
        adderror('sdf:syntax:incomplete', ...
                 'The field ''variables'' (lowercase) is not specified'); 
    end
    if ~ismember('factors', fields)
        if ismember('transform', fields) || model.options.ImplicitFactors
            model.factors = {};
        else 
            adderror('sdf:syntax:incomplete', ...
                     'The field ''factors'' (lowercase) is not specified'); 
        end
    end
    if ~ismember(fields, 'factorizations')
        adderror('sdf:syntax:incomplete', ...
                 'The field ''factorizations'' (lowercase) is not specified'); 
    end

    %% Add book keeping fields if necessary
    if ~isfield(model, 'names')
        model.names = struct;
    end 

    %% Check variable structure
    if ~iscell(model.variables) && ~isstruct(model.variables)
        model.variables = {model.variables};
    end
    if ~iscell(model.variables) % hence a struct
        model.names.variables = fieldnames(model.variables);
        model.variables = struct2cell(model.variables);
    end 
    if isfield(model.names, 'variables')
        ia = cellfun(@(s) any(cellfun(@(t) strcmp(s,t), reservedKeywords)), ...
                     model.names.variables);
        ia = find(ia);
        if ~isempty(ia)
            for k = 1:length(ia)
                adderror('variables', ia(k), 'sdf:syntax:reservedKeyword', ...
                         ['The name ''%s'' is reserved and cannot be used as variable ' ...
                          'name'], ...
                         model.names.variables{ia(k)});
            end
        end
    end
    model.variables = model.variables(:).';

    %% Check variable properties
    for k = 1:length(model.variables)
        % variables should be (potentially nested) numerical values
        if ~checknumeric(model.variables{k})
            adderror('variables', k, 'sdf:syntax:notNumeric', ...
                     ['The variable %s is not numerical, or a potentially nested ' ...
                      'cell of numerical values.'], ...
                     getvariablename(k));
        end
    end
   
    %% Check factor structure
    % track some syntax information
    factorstructure.isarray = false;
    
    % if not a cell, wrap
    if ~iscell(model.factors) && ~isstruct(model.factors)
        % auto wrapping 
        if isnumeric(model.factors)
            model.factors = num2cell(model.factors);
        else 
            model.factors = {model.factors};
        end
        factorstructure.isarray = true;
    end
    if ~iscell(model.factors) % hence a struct
        model.names.factors = fieldnames(model.factors);
        model.factors = struct2cell(model.factors);
    end 
    model.factors = model.factors(:).';

    % Check names
    if isfield(model.names, 'factors') && model.options.ImplicitFactors
        IFnames = ~cellfun(@isempty, regexp(model.names.factors, '^IF[0-9]+$'));
        if any(IFnames)
            IFnames = find(IFnames, 1);
            adderror('factors', IFnames, 'sdf:syntax:reservedName', ...
                     ['If implicit factors are allowed, factor names cannot have ' ...
                      'the form ''IFx'' in which x is a number.']);
        end
    end
    
    % Wrap if a single factor with functions
    funcs = cellfun(isfunc, model.factors);
    if length(model.factors) > 1 && all(funcs(2:end)) && ...
            ((isscalar(model.factors{1}) && isnumeric(model.factors{1})) || ...
             ischar(model.factors{1}))
        model.factors = {model.factors};
    end

    for k = 1:length(model.factors)
        % case 1: non-transformed variable => wrap
        if ~iscell(model.factors{k})
            model.factors{k} = model.factors(k);
        end
        % case 2: single subfactor => wrap
        funcs = cellfun(isfunc, model.factors{k});
        if any(funcs)
            if sum(funcs) ~= length(funcs) - 1 
                % all but the first entries should be function handles, except
                % when there is an explicit conversion to a constant
                if sum(funcs) == length(funcs) - 2 && ...
                        strcmp(model.factors{k}{2}, '(constant)')
                    model.factors{k} = model.factors(k);
                elseif ~isnumeric(model.factors{k}{1}) && ~ ...
                        ischar(model.factors{k}{1})
                    adderror('factors', k, 'sdf:syntax:subfactorConversion', ...
                             ['Could not convert factor %s into subfactors: ' ...
                              'the first entry should be a variable reference.'], ...
                             getfactorname(k));
                else 
                    adderror('factors', k, 'sdf:syntax:subfactorConversion', ...
                             ['Could not convert factor %s into subfactors: a variable ' ...
                              'can only be followed by transformations'], ...
                             getfactorname(k));
                end
            elseif iscell(model.factors{k}{1})
                % first entry cannot be a cell (may change in future)
                adderror('factors', k, 'sdf:syntax:subfactorConversion', ...
                         ['Could not convert factor %s into subfactors: the input ' ...
                          'of a transformation should be a variable'], ...
                         getfactorname(k));
            else % valid => wrap
                model.factors{k} = model.factors(k);
            end
        end
        % case 3a: group explicit constants
        strs = find(cellfun(@ischar, model.factors{k}));
        constants = strs(cellfun(@(s) strcmp(s, '(constant)'), ...
                                 model.factors{k}(strs)));
        for l = constants(:).'
            % merge constants
            model.factors{k}{l-1} = {model.factors{k}{l-1} model.factors{k}{l}};
        end
        % remove constants
        model.factors{k}(constants) = [];

        % case 3b: subfactors but not a cell => wrap
        noncellidx = find(~cellfun(@iscell, model.factors{k}));
        for l = noncellidx(:).'
            model.factors{k}{l} = model.factors{k}(l);
        end
    end

    %% Dereference variable references, find constants and check validity names
    if ~isfield(model, 'isconstant')
        model.isconstant = cell(1, length(model.factors));
    end
    for k = 1:length(model.factors) % factors
        model.isconstant{k} = false(size(model.factors{k}));
        for l = 1:numel(model.factors{k}) % subfactors
            if length(model.factors{k}{l}) > 1 && ischar(model.factors{k}{l}{2}) ...
                    && strcmp(model.factors{k}{l}{2}, '(constant)');
                % explicitly defined constant
                model.isconstant{k}(l) = true;
                model.factors{k}{l}(2) = [];
            elseif ischar(model.factors{k}{l}{1}) % named reference
                if ~isfield(model.names, 'variables')
                    adderror('factors', k, 'sdf:syntax:unknownVariable', ...
                             'Factor %s depends on the unknown variable %s', ...
                             getfactorname(k), model.factors{k}{l}{1});
                end
                idx = findname(model.factors{k}{l}{1}, model.names.variables);
                if isempty(idx)
                    adderror('factors', k, 'sdf:syntax:unknownVariable', ...
                             'Factor %s depends on the unknown variable %s', ...
                             getfactorname(k), model.factors{k}{l}{1});
                else
                    model.factors{k}{l}{1} = idx;
                end
            elseif isscalar(model.factors{k}{l}{1}) && isnumeric(model.factors{k}{l}{1}) % indexed reference
                idx = model.factors{k}{l}{1};
                if ~(round(idx) == idx && idx > 0 && idx <= ...
                     length(model.variables))
                    adderror('factors', k, 'sdf:syntax:unknownVariable', ...
                             'Factor %s depends on the unknown variable %d', ...
                             getfactorname(k), model.factors{k}{l}{1});
                end
            else % implicitly defined constant
                model.isconstant{k}(l) = true;
            end
        end
    end

    %% Check constants 
    for k = 1:length(model.factors)
        for l = find(model.isconstant{k}(:).')
            % variables should be (potentially nested) numerical values
            if ~checknumeric(model.factors{k}{l}{1})
                adderror('factors', k, 'sdf:syntax:notNumeric', ...
                         ['A constant in factor %s is not numerical, or a potentially ' ...
                          'nested cell of numerical values.'], ...
                         getfactorname(k));
            end        
        end
    end
    
    % If transformations are disabled, check if no transformations are
    % defined
    if options.NoTransformations
        for k = 1:length(model.factors)
            if any(cellfun(@(f) length(f) > 1, model.factors{k}))
                adderror('factors', k, 'sdf:syntax:disabledFeatureUsed', ...
                         ['No transformations can be used when the option ' ...
                          '''NoTransformations'' is given.']);
            end
        end
    end
    if options.OnlySimpleFactors
        idx = cellfun(@(f) length(f) == 1, model.factors);
        if ~all(idx)
            adderror('factors', find(idx,1), 'sdf:syntax:disabledFeatureUsed', ...
                     ['Variables cannot be concatenated when the option ' ...
                      '''OnlySimpleFactors'' is given.']);
        end
    end
        
    %% Check transform structure
    % structure: { {varidx}, {factoridx}, func1, arg1, arg2, func2, ...}
    % if factoridx is omitted, varidx is taken
    
    % Book keeping
    isemptyfactor = @(f) length(f) == 1 && length(f{1}) == 1 && ...
        isempty(f{1}{1});
    if isfield(model, 'factors')
        % dectect empty factors, as they can be overwritten
        emptyfactors = cellfun(isemptyfactor, model.factors);
        cache.factors.origin = zeros(1, length(emptyfactors));
        cache.factors.origin(~emptyfactors) = -1;
    else 
        emptyfactors = [];
        cache.factors.origin = [];
    end
    
    if isfield(model, 'transform') && ~isempty(model.transform)
        if isstruct(model.transform)
            model.names.transform = fieldnames(model.transform);
            model.transform = struct2cell(model.transform);
            if ~isempty(model.transform)
                for k = find(cellfun(@(t) ~iscell(t), model.transform(:)'))
                    if isnumeric(model.transform{k}) || ...
                            ischar(model.transform{k})
                        % auto wrap
                        model.transform{k} = model.transform(k);
                    else 
                        adderror('transform', k, 'sdf:syntax:invalidModel', ...
                                 ['A transform entry should be a cell, a ' ...
                                  'numerical array, or a string.']);
                    end
                end
            end
        end
        % if not a cell, doubel wrapping can maybe be applied
        if ~iscell(model.transform)
            if isnumeric(model.transform) || ischar(model.transform)
                model.transform = {{model.transform}};
            else 
                adderror('transform', idx, 'sdf:syntax:invalidModel', ...
                         'The transform field should be a struct or cell of cells');
            end
        end
        if any(~cellfun(@iscell, model.transform))
            % auto-wrap if unambiguous, e.g., if a function handle is present
            if any(cellfun(@(f) isa(f, 'function_handle'), model.transform))
                model.transform = {model.transform(:).'};
            else 
                idx = find(~cellfun(@iscell, model.transform(:)'),1);
                adderror('transform', idx, 'sdf:syntax:invalidModel', ...
                         'The transform field should be a struct or cell of cells');
            end 
        end
        
        model.transform = model.transform(:)';
        nbfactors = length(emptyfactors);

        for k = 1:length(model.transform)
            % each entry should be a cell
            if ~iscell(model.transform{k})
                adderror('transform', k, 'sdf:syntax:invalidModel', ...
                         'Each transform entry should be a cell');
            end
            model.transform{k} = model.transform{k}(:).';
            % If the first entry (variables) is not a cell, auto-wrap
            if ~iscell(model.transform{k}{1}) 
                if isnumeric(model.transform{k}{1})
                    model.transform{k}{1} = num2cell(model.transform{k}{1}(:).');
                elseif ischar(model.transform{k}{1})
                    model.transform{k}{1} = model.transform{k}(1);
                else % should never happen
                    adderror('transform', {k, 1}, 'sdf:syntax:unknownVariable', ...
                             ['The first entry of a transform operation should ' ...
                              'be a cell of valid variable references (string ' ...
                              'names, or indices).']);
                end
            end
            model.transform{k}{1} = model.transform{k}{1}(:).';
            % if no factor references are given, copy variable references
            if length(model.transform{k}) == 1 || ...
                    isfunc(model.transform{k}{2})
                model.transform{k} = model.transform{k}([1 1 2:end]);
                fvidx = 1; % factor variable index needed for error messages
            else
                fvidx = 2;
            end
            % auto-wrap factor references
            if ~iscell(model.transform{k}{2}) 
                if ischar(model.transform{k}{2})
                    model.transform{k}{2} = model.transform{k}(2);
                elseif isnumeric(model.transform{k}{2})
                    model.transform{k}{2} = num2cell(model.transform{k}{2});
                else  % should never happend
                    adderror('transform', {k, fvidx}, 'sdf:syntax:invalidReference', ...
                             'Factor references should be either strings or integers');
                end
            end

            % step 1: Resolve variable references
            strrefs = cellfun(@ischar, model.transform{k}{1});
            intrefs = cellfun(@(r) isscalar(r) && isnumeric(r) , model.transform{k}{1});
            if ~all(strrefs | intrefs)
                adderror('transform', {k, 1}, 'sdf:syntax:invalidReference', ...
                         ['A variable reference should be either a string or ' ...
                          'an integer.']);
            end
            if any(strrefs) && ~isfield(model.names, 'variables')
                adderror('transform', {k, 1}, 'sdf:syntax:invalidReference', ...
                         ['A variable reference can only be a string if ' ...
                          'model.variables is a struct.']);
            end
            for l = find(strrefs)
                idx = findname(model.transform{k}{1}{l}, ...
                               model.names.variables);
                if isempty(idx)
                    adderror('transform', {k, 1}, 'sdf:syntax:unknownVariable', ...
                             'No variable named %s is known.', ...
                             model.transform{k}{1}{l});
                else
                    model.transform{k}{1}{l} = idx;
                end
            end
            % check integer references
            for l = find(intrefs)
                idx = model.transform{k}{1}{l};
                if ~(round(idx) == idx && idx > 0 && idx <= ...
                     length(model.variables))
                    adderror('transform', {k, 1}, 'sdf:syntax:unknownVariable', ...
                             'No variable with id %d is known.', ...
                             model.transform{k}{1}{l});
                end
            end

            % step 2: check factor references
            % lengths should match
            if length(model.transform{k}{1}) ~= ...
                        length(model.transform{k}{2})
                adderror('transform', {k, 1:2}, 'sdf:syntax:dimensionMismatch', ...
                         ['The length of the variable references should ' ...
                          'match the length of the factor references']);
            end
            
            % either all strings or all integers
            allstrings = all(cellfun(@ischar, model.transform{k}{2}));
            allints    = all(cellfun(@isnumeric, model.transform{k}{2}));
            if isempty(model.factors) && allstrings
                model.names.factors = {};
            end
            if isfield(model.names, 'factors') && ~allstrings
                adderror('transform', {k, fvidx}, 'sdf:syntax:invalidReference', ...
                         ['Only string references can be used if named factors ' ...
                          'are used (model.factors is a struct).']);
            elseif ~isfield(model.names, 'factors') && ~allints
                adderror('transform', {k, fvidx}, 'sdf:syntax:invalidReference', ...
                         ['Only integer references can be used if no named factors ' ...
                          'are used (model.factors is a cell.']);
            end
            % references made to defined factors
            if allstrings
                [~, ia] = intersect(model.transform{k}{2}, ...
                                    model.names.factors);
                if ~isempty(ia)
                    for l = 1:length(ia)
                        if ia(l) <= nbfactors
                            adderror('transform', {k,fvidx}, ...
                                     'sdf:syntax:overwriteFactor', ...
                                     ['Transform %s cannot overwrite ' ...
                                      'factor %s.'], ...
                                     gettransformname(k), ...
                                     getfactorname(ia(l)));
                        else 
                            adderror('transform', {k,fvidx}, ...
                                     'sdf:syntax:overwriteFactor', ...
                                     ['Transform %s cannot overwrite ' ...
                                      'the factor %s defined by a previous ' ...
                                      'transform.'], ...
                                     gettransformname(k), ...
                                     getfactorname(ia(l)));
                        end
                    end
                end
                % if no errors, convert strings to numbers
                [factornames,~,ic] = unique(model.transform{k}{2}, 'stable');
                if length(factornames) == length(model.transform{k}{2})
                    ic = ic + length(model.names.factors);
                    model.names.factors = [model.names.factors; factornames(:)];
                    model.transform{k}{2} = num2cell(ic);
                else 
                    adderror('transform', {k,fvidx}, 'sdf:syntax:overwriteFactor', ...
                             ['Transform %s defines some factors multiple ' ...
                              'times.'], ...
                              gettransformname(k));
                end
            else 
                ia = cell2mat(model.transform{k}{2}) <= length(model.factors);
                ia = [model.transform{k}{2}{ia}];
                for l = ia
                    if l <= nbfactors && ~emptyfactors(l)
                        % overwrites a defined factor
                        adderror('transform', {k,fvidx}, ...
                                 'sdf:syntax:overwriteFactor', ...
                                 ['Transform %s cannot overwrite factor ' ...
                                  '%d.'], ...
                                 gettransformname(k), l);
                    elseif ~(isempty(model.factors{l}) || ...
                             isemptyfactor(model.factors{l}))
                        % overwrites result of previous transform
                        adderror('transform', {k,fvidx}, ...
                                 'sdf:syntax:overwriteFactor', ...
                                 ['Transform %s cannot overwrite ' ...
                                  'the factor %d defined by a previous ' ...
                                  'transform.'], ...
                                 gettransformname(k), l);
                    end
                end
                if any(diff(sort(cell2mat(model.transform{k}{2})))==0)
                    adderror('transform', {k,fvidx}, 'sdf:syntax:overwriteFactor', ...
                             ['Transform %s defines some factors multiple ' ...
                              'times.'], ...
                             gettransformname(k));
                end
            end
            % no errors means no factors are overwritten    
            
            if options.NoTransformations && length(model.transform{k}) >= 3
                adderror('transform', k, 'sdf:syntax:disabledFeatureUsed', ...
                         ['No transformations can be used when the options ' ...
                          '''NoTransformations'' is given.']);
            end

            % step 3: group functions
            if length(model.transform{k}) >= 3 && ~ ...
                    isfunc(model.transform{k}{3})
                adderror('transform', {k, fvidx+1}, ...
                         'sdf:syntax:invalidTransformation', ...
                         'Transform %s contains an invalid transformation', ...
                         gettransformname(k));
            end
            
            funcidx = find(cellfun(isfunc, model.transform{k}));
            nbfuncs = length(funcidx);
            funcs = cell(length(model.transform{k}{2}), nbfuncs);
            
            for l = 1:nbfuncs
                f = model.transform{k}{funcidx(l)};
                if l < nbfuncs
                    args = model.transform{k}(funcidx(l)+1:funcidx(l+1)-1);
                else
                    args = model.transform{k}(funcidx(l)+1:end);
                end
                % convert arguments
                for m = 1:length(args)
                    if iscell(args{m})
                        if numel(args{m}) == 1
                            args{m} = repmat(args{m}, length(model.transform{k}{2}), 1);
                        elseif numel(args{m}) == ...
                                length(model.transform{k}{2})
                            args{m} = args{m}(:);
                        else 
                            adderror('transform', {k, funcidx(l)+m+1-fvidx}, ...
                                     'sdf:syntax:invalidParameter', ...
                                     ['Cell type parameters should have either ' ...
                                      'length one or the number of factors to ' ...
                                      'transform.']);
                        end
                    else 
                        args{m} = repmat(args(m), length(model.transform{k}{2}), 1);
                    end
                end
                % create functions
                for m = 1:length(model.transform{k}{2})
                    % ugly hack to display right function name
                    arg = cellfun(@(a) a{m}, args, 'UniformOutput', false); %#ok
                    eval(sprintf('funcs{m,l} = @(z,task) %s(z, task, arg{:});', ...
                                 func2str(f)));
                end
            end
            
            % step 4: create factors
            for l = 1:length(model.transform{k}{2})
                vidx = model.transform{k}{1}{l};
                idx  = model.transform{k}{2}{l};
                model.factors{idx} = {{vidx, funcs{l,:}}};
            end
            
            % step 5: keep track of which transformation added a factor
            cache.factors.origin(cat(2, model.transform{k}{2}{:})) = k;
        end

        % Fix previously constant factors that have been transformed
        constants = cellfun(@(f) numel(f) == 1 && f == 1, model.isconstant);
        for k = find(constants)
            if emptyfactors(k) % was empty before
                model.isconstant{k} = 0;
            end
        end
        % add non-constant for new fields
        for k = length(model.isconstant)+1:length(model.factors);
            model.isconstant{k} = false;
        end
    end
        
    %% Check factorization structure
    factorizationtype = 'cell';
    if ~isstruct(model.factorizations) && ~iscell(model.factorizations)
        adderror('factorizations', k, 'sdf:syntax:invalidModel', ...
                 'model.factorizations should be either a struct or a cell');
    end
    if ~iscell(model.factorizations)
        if numel(model.factorizations) == 1 % no struct array
            model.names.factorizations = fieldnames(model.factorizations);
            model.factorizations = struct2cell(model.factorizations);
            factorizationtype = 'struct';
        else 
            fn = fieldnames(model.factorizations);
            vals = squeeze(struct2cell(model.factorizations));
            model.factorizations = cell(1, numel(model.factorizations));
            for k = 1:numel(model.factorizations)
                tmp = [fn(:).'; vals(:,k).'];
                model.factorizations{k} = struct(tmp{:});
            end
            factorizationtype = 'structarray';
        end
    end
    model.factorizations = model.factorizations(:).';

    for k = 1:length(model.factorizations)
        % each factorization should be a struct
        if ~isstruct(model.factorizations{k})
            adderror('factorizations', k, 'sdf:syntax:invalidModel', ...
                     'Each factorization should be a struct.')
        end
        if isfield(model.factorizations{k}, 'type') || ...
                isfield(model.factorizations{k}, 'factors')
            if ~isfield(model.factorizations{k}, 'type')
                adderror('factorizations', k, 'sdf:syntax:unknownFactorization', ...
                         'No factorization found for factorization %s', ...
                         getfactorizationname(k));
            end
            if ~isfield(model.factorizations{k}, 'factors')
                adderror('factorizations', k, 'sdf:syntax:noFactors', ...
                         'No factors found for factorization %s', ...
                         getfactorizationname(k));
            end

            % check if valid type
            if ~ischar(model.factorizations{k}.type)
                 adderror('factorizations', {k,'type'}, 'sdf:syntax:invalidModel', ...
                         'The factorization type of factorization %s should be a string', ...
                          getfactorizationname(k));               
            end
            if ~any(strcmp(model.factorizations{k}.type, ...
                           knownfactorizations))
                adderror('factorizations', {k,'type'}, 'sdf:syntax:unknownFactorization', ...
                         'No known factorization type found for factorization %s', ...
                         getfactorizationname(k));
            end 
            % Check if valid factors
            model.factorizations{k}.factors = ...
                fmtfactors(model.factorizations{k}.factors, model.factorizations{k}.type);    
            if ~iscell(model.factorizations{k}.factors)
                 adderror('factorizations', {k,'factors'}, 'sdf:syntax:incorrectStructure', ...
                         'The factors field of factorization %s should be a cell of factor references', ...
                          getfactorizationname(k));               
            end
            % check if only type
            fn = fieldnames(model.factorizations{k});
            others = ~strcmp('data',fn) & ~strcmp('type',fn) & ...
                          ~strcmp('factors', fn);
            fn = fn(others);
            type = cellfun(@(s) any(strcmp(s, knownfactorizations)), fn);
            type = fn(type);
            
            if ~isempty(type)
                if length(type) > 1 || ~strcmp(type, model.factorizations{k}.type)
                    adderror('factorizations', k, 'sdf:syntax:multipleFactorizations', ...
                             'Factorization %s can only have one factorization type', ...
                             getfactorizationname(k));
                end
                % if type is equal, check if all references are the same
                % first format the references
                f1 = fmtfactors(model.factorizations{k}.factors, type{1});
                f2 = fmtfactors(model.factorizations{k}.(type{1}), type{1});
                % then derefence
                f1 = dereffacts(k, f1);
                f2 = dereffacts(k, f2);

                f1 = cat(2,f1{:});
                f2 = cat(2,f2{:});
                if any(strcmpi(type{1}, {'btd', 'lmlra'}))
                    f1 = cat(2,f1{:});
                    f2 = cat(2,f2{:});
                end
                % now check
                if (length(f1) ~= length(f2) || any(f1(:) ~= f2(:)))
                    adderror('factorizations', {k,'factors'}, 'sdf:syntax:multipleFactorizations', ...
                             ['Factorization %s defines two %s factorizations with conflicting factors ' ...
                              'in the factors field and the %s field.'], ...
                             getfactorizationname(k),type{1},type{1}); 
                end
                % remove unnecessary field
                tmp = struct2cell(model.factorizations{k});
                fields = fieldnames(model.factorizations{k});
                mask = ~cellfun(@(s) strcmp(s, type{1}), fields);
                model.factorizations{k} = cell2struct(tmp(mask), fields(mask));                                
            end        
        else 
            fn = fieldnames(model.factorizations{k});
            % check type
            model.factorizations{k}.type = find(~strcmp('data',fn));
            if isempty(model.factorizations{k}.type)
                adderror('factorizations', k, 'sdf:syntax:unknownFactorization', ...
                         'No factorization found for factorization %s', ...
                         getfactorizationname(k));
            else
                type = intersect(fn(model.factorizations{k}.type), ...
                                 knownfactorizations);
                if length(type) == 1
                    model.factorizations{k}.type = type{1};
                elseif isempty(type)
                    % change: unknown
                    adderror('factorizations', k, 'sdf:syntax:unknownFactorization', ...
                             'No known factorization type found for factorization %s', ...
                             getfactorizationname(k));
                elseif length(type) > 1
                    adderror('factorizations', k, 'sdf:syntax:multipleFactorizations', ...
                             'Factorization %s can only have one factorization type', ...
                             getfactorizationname(k));            
                end

            end
            
            % extract factors
            model.factorizations{k}.factors = ...
                model.factorizations{k}.(model.factorizations{k}.type);
            type = model.factorizations{k}.type;
            % remove the type field (faster than rmfield)
            tmp = struct2cell(model.factorizations{k});
            fields = fieldnames(model.factorizations{k});
            mask = ~cellfun(@(s) strcmp(s, type), fields);
            model.factorizations{k} = cell2struct(tmp(mask), fields(mask));
        end
        
        % Is the factorization allowed?
        type = model.factorizations{k}.type;
        if ~any(cellfun(@(s) strcmp(type,s), allowedfactorizations))
            adderror('factorizations', k, 'sdf:syntax:disabledFeatureUsed', ...
                     'The use of the factorization ''%s'' is disabled', ...
                     type);
        end
        
        % reformat factors if necessary
        model.factorizations{k}.factors = fmtfactors(model.factorizations{k}.factors, ...
                                                     model.factorizations{k}.type);
    end
    
    %% Check factor references format
    isref = @(r) (isscalar(r) && isnumeric(r)) || ischar(r);
    for k = 1:length(model.factorizations)
        factors = model.factorizations{k}.factors;
        location = {k, model.factorizations{k}.type};
        switch model.factorizations{k}.type
          case {'cpd', 'cpdi', 'll1', 'regL1', 'regL2'}
            % correct format is a cell of integer or string references
            if ~iscell(factors)
                adderror('factorizations', location, 'sdf:syntax:incorrectStructure', ...
                         'The factors field for factorization %s is not a cell', ...
                         getfactorizationname(k));
            elseif ~all(cellfun(isref, factors))
                if strcmp(model.factorizations{k}.type, 'll1') && ... 
                        all(cellfun(@iscell, factors)) && ...
                        all(cellfun(@(f) all(cellfun(isref,f)), factors))
                    adderror('factorizations', location, 'sdf:syntax:wrongTypeFactors', ...
                             ['In factorization %s, btd type factors are ' ...
                              'detected, while ll1 currently only accepts CPD ' ...
                              'type factors combined with an L field.'], ...
                             getfactorizationname(k));  
                end              
                adderror('factorizations', location, 'sdf:syntax:incorrectStructure', ...
                         ['One or more entries in the factors field cell of ' ...
                          'factorization %s is not a valid factor reference'], ...
                         getfactorizationname(k));            
            end
          case {'btd', 'lmlra'}
            % correct format is a cell of cells of integer references
            if ~iscell(factors) || any(~cellfun(@iscell, factors))
                adderror('factorizations', location, 'sdf:syntax:incorrectStructure', ...
                         ['The factors field for factorization %s is not a cell ' ...
                          'of cells'], ...
                         getfactorizationname(k));
            elseif ~all(cellfun(@(c) all(cellfun(isref, c)), factors))
                adderror('factorizations', location, 'sdf:syntax:incorrectStructure', ...
                         ['One or more entries in the factors field cell of ' ...
                          'factorization %s is not a valid factor reference'], ...
                         getfactorizationname(k));            
            end        
        end
    end
    
    %% Dereference factor references
    for k = 1:length(model.factorizations)
        model.factorizations{k}.factors = ... 
            dereffacts(k, model.factorizations{k}.factors);
    end

    %% Check data fields
    for k = 1:length(model.factorizations)
        factors = model.factorizations{k}.factors;
        switch model.factorizations{k}.type
          case {'cpd', 'cpdi', 'll1', 'btd', 'lmlra'}
            if ~isfield(model.factorizations{k}, 'data')
                adderror('factorizations', k, 'sdf:syntax:noDataFound', ...
                         'Factorization %s has no data field as required for %s', ...
                         getfactorizationname(k), ...
                         model.factorizations{k}.type);
            end
          case {'regL0', 'regL1', 'regL2'}
            if isfield(model.factorizations{k}, 'data')
                data = model.factorizations{k}.data;
                location = {k, 'data'};
                if ~iscell(data) % auto wrap
                    model.factorizations{k}.data = {data};
                end 
            end
        end
    end

    %% Check extra fields
    for k = 1:length(model.factorizations)
        switch model.factorizations{k}.type
          case 'll1'
            % ll1 should have an L field
            if ~isfield(model.factorizations{k}, 'L')
                adderror('factorizations', k, 'sdf:syntax:noL', ...
                         'Factorization %s does not specify the field L', ...
                         getfactorizationname(k));      
            end
            % and it should be a vector
            if ~isvector(model.factorizations{k}.L)
                adderror('factorizations', {k, 'L'}, 'sdf:syntax:incorrectL', ...
                         'Factorization %s: L should be a vector', ...
                         getfactorizationname(k));      
            end
            model.factorizations{k}.L = model.factorizations{k}.L(:).';
        end

        % check additional fields
        if isfield(model.factorizations{k}, 'issymmetric')
            if ~isscalar(model.factorizations{k}.issymmetric) || ...
                    ~islogical(model.factorizations{k}.issymmetric)
                adderror('factorizations', {k,'issymmetric'}, ...
                         'sdf:syntax:invalidValue',...
                         ['Factorization %s: the issymmetric option should ' ...
                          'have true or false as value.'], ...
                         getfactorizationname(k));
            end
        end
        if isfield(model.factorizations{k}, 'weight')
            if ~isscalar(model.factorizations{k}.weight) || ...
                    ~isnumeric(model.factorizations{k}.weight) || ...
                    isnan(model.factorizations{k}.weight)
                adderror('factorizations', {k,'weight'}, ...
                         'sdf:syntax:invalidValue',...
                         ['Factorization %s: the weight options should be a ' ...
                          'scalar value.'], getfactorizationname(k));
            end
        end
        if isfield(model.factorizations{k}, 'relweight')
            if ~isscalar(model.factorizations{k}.relweight) || ... 
                    ~isnumeric(model.factorizations{k}.relweight) || ...
                    isnan(model.factorizations{k}.relweight)
                adderror('factorizations', {k,'relweight'}, ...
                         'sdf:syntax:invalidValue',...
                         ['Factorization %s: the relweight options should be a ' ...
                          'scalar value.'], getfactorizationname(k));
            end
        end
    end

    % Add missing (rel)weights if necessary
    hasweights    = cellfun(@(f) isfield(f, 'weight'),    model.factorizations);
    hasrelweights = cellfun(@(f) isfield(f, 'relweight'), model.factorizations);
    if any(hasweights) && any(hasrelweights)
        adderror('factorizations', [], 'sdf:syntax:mixedWeights', ...
                 ['All factorizations should define the same weight type, i.e., ' ...
                  'either all factorizations have a weight, or all factorizations ' ...
                  'have a relweight.']);
    elseif any(hasweights)
        idx = find(~hasweights);
        for k = idx(:).'
            model.factorizations{k}.weight = 1;
        end
    elseif any(hasrelweights)
        idx = find(~hasrelweights);
        for k = idx(:).'
            model.factorizations{k}.relweight = 1;
        end
    end
    
    %% Check if consistent
    if ~isconsistent
        varargout = {model, errorlist};
        return;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 2: Consistency
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.Print, 
        options.StopOnError = false; 
        options.PrintErrors = false;
    end

    %% Expand and check all factors
    cache.factors.expanded = cell(1, length(model.factors));
    cache.factors.sequance = cell(1, length(model.factors));
    cache.factors.errors = cell(1, length(model.factors));
    cache.factors.suberrors = cell(1, length(model.factors));
    for k = 1:length(model.factors)
        cache.factors.sequence{k} = cell(size(model.factors{k}));
        cache.factors.suberrors{k} = cell(size(model.factors{k}));
        % determine the origin: factor or transform?
        if cache.factors.origin(k) == -1
            origin = 'factors';
            originid = k;
            originmsg = sprintf('factor %s', getfactorname(k));
        elseif cache.factors.origin(k) > 0 
            origin = 'transform';
            originid = cache.factors.origin(k);
            originmsg = sprintf('factor %s (generated by transform %s)', ...
                                getfactorname(k), gettransformname(originid));
        else 
            origin = 'factor';
            originid = nan;
            originmsg = sprintf('factor %s', getfactorname(k));
        end
        for l = 1:numel(model.factors{k})
            cache.factors.sequence{k}{l} = cell(1, length(model.factors{k}{l}));
            % get variables (should not give errors for correct syntax)
            if model.isconstant{k}(l) 
                cache.factors.sequence{k}{l}{1} = model.factors{k}{l}{1};
            else
                cache.factors.sequence{k}{l}{1} = ...
                    model.variables{model.factors{k}{l}{1}};
            end
            % apply transformations from left to right
            cache.factors.suberrors{k}{l} = cell(1, length(model.factors{k}{l})-1);
            for m = 2:length(model.factors{k}{l})
                try 
                    cache.factors.sequence{k}{l}{m} = ...
                        model.factors{k}{l}{m}(cache.factors.sequence{k}{l}{m-1},[]);
                catch e
                    %FIXME
                    % if strcmpi(origin, 'factors')
                    %     if nbfactors == 2
                    %         tmpid = {origid};
                    %     end
                    % end 
                    inchecker = strcmp('sdf_check', e.stack(1).name) || ...
                        strncmp('sdf_check', e.stack(2).name, 9);
                    switch e.identifier
                      case 'MATLAB:UndefinedFunction'
                        if inchecker
                            adderror(origin, originid, 'sdf:consistency:unknownTransformation', ...
                                     'In %s, the transformation %s is unknown', ...
                                     originmsg, ...
                                     getfunname(model.factors{k}{l}{m}));
                        else 
                            adderror(origin, originid, 'sdf:consistency:structError', ...
                                     ['While evaluating %s, the transformation ' ...
                                      '%s threw the %s error:\n\n    %s'], ...
                                     originmsg, ...
                                     getfunname(model.factors{k}{l}{m}), ...
                                     e.identifier, e.message);
                        end
                        cache.factors.suberrors{k}{l}{m-1} = ...
                            errorlist{end};
                      otherwise
                        adderror(origin, originid, 'sdf:consistency:structError', ...
                                 ['While evaluating %s, the transformation ' ...
                                  '%s threw the %s error:\n\n    %s'], ...
                                 originmsg, ...
                                 getfunname(model.factors{k}{l}{m}), ...
                                 e.identifier, e.message);
                        cache.factors.suberrors{k}{l}{m-1} = errorlist{end};
                    end
                end
            end 
        end
        % concatenate results if necessary
        if numel(model.factors{k}) == 1
            cache.factors.expanded{k} = cache.factors.sequence{k}{1}{end};
        else 
            tmp = cellfun(@(f) f{end}, cache.factors.sequence{k}, 'UniformOutput', ...
                          false);
            try 
                cache.factors.expanded{k} = cell2mat(tmp);
            catch e
                if strcmp(e.identifier, 'MATLAB:catenate:dimensionMismatch')
                    adderror('factors', k, 'sdf:consistency:concatError', ...
                             ['Could not concatenate subfactors in Factor %s: the ' ...
                              'dimensions do not match'], ...
                             getfactorname(k));
                    cache.factors.errors{k}{end+1} = errorlist{end};
                else 
                    rethrow(e);
                end
            end
        end 
    end

    %% Check factorizations
    cache.factorizations.errors = cell(1, length(model.factorizations));
    for k = 1:length(model.factorizations)
        % Step 1: Check data
        location = {k, 'data'};
        switch model.factorizations{k}.type
          case {'cpd', 'cpdi', 'll1', 'btd', 'lmlra'}
            data = model.factorizations{k}.data;
            % is it a known type?
            try 
                type = getstructure(data);
            catch e
                adderror('factorizations', location, 'sdf:consistency:unknownDataType', ...
                         ['The data field of factorization %s contains an ' ...
                          'unknown data type'], ...
                         getfactorizationname(k));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
                continue;
            end
            % in the case of cpdi, it should be an incomplete tensor
            if strcmpi(model.factorizations{k}.type, 'cpdi')
                if strcmpi(type, 'full')
                    % try to format if full
                    data = fmt(data);
                    if ~isstruct(data) || ~data.incomplete
                        adderror('factorizations', location, 'sdf:consistency:invalidDataType',...
                                 ['The ''cpdi'' factorization only accepts ' ...
                                  'incomplete data.']);
                    end
                    type = 'incomplete';
                    model.factorizations{k}.data = data;
                elseif ~strcmpi(type, 'incomplete')
                    adderror('factorizations', location, 'sdf:consistency:invalidDataType',...
                             ['The ''cpdi'' factorization only accepts ' ...
                              'incomplete data.']);
                end 
            end
            % if yes, is the data well specified?
            if ~isvalidtensor(data)
                cmd = sprintf('isvalidtensor(%s, true)', ...
                              getsyntaxline('factorizations', k, 'data'));
                adderror('factorizations', location, 'sdf:consistency:incorrectDataFormat', ...
                         ['The data of factorization %s is has an incorrect ' ...
                          'structure. Use %s to check for errors.'], ...
                         getfactorizationname(k), createlink(cmd));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
                continue;
            end
            model.factorizations{k}.datatype = type;
            unstructuredtypes = {'full', 'incomplete', 'sparse'};
            model.factorizations{k}.isstructured = ~any(cellfun(@(s) strcmpi(type, ...
                                                              s), unstructuredtypes));          
          case {'regL0', 'regL1', 'regL2'}
            if isfield(model.factorizations{k}, 'data')
                data = model.factorizations{k}.data;
                if any(~cellfun(@isnumeric, data))
                    adderror('factorizations', location, 'sdf:consistency:incorrectDataFormat', ...
                             ['Factorization %s contains one or more non-numerical ' ...
                              'data fields'], ...
                             getfactorizationname(k));
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                    continue;
                elseif length(factors) ~= length(data)
                    adderror('factorizations', location, 'sdf:consistency:incorrectNumberDatasets', ...
                             ['The number of datasets (%d) differs from the number ' ...
                              'of factors (%d) in factorization %s'], ...
                             length(data), length(factors), ...
                             getfactorizationname(k));
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                    continue;
                end
            end
        end        
        
        % Step 2: check factors
        location = {k, model.factorizations{k}.type};
        switch model.factorizations{k}.type
          case {'cpd', 'cpdi'}
            U = cache.factors.expanded([model.factorizations{k}.factors{:}]);
            % factor matrices should be matrices
            if any(~cellfun(@ismatrix, U)) || any(~cellfun(@ isnumeric, U))
                adderror('factorizations', location, 'sdf:consistency:wrongTypeFactors', ...
                         'Factorization %s: in a cpd all factors should be numerical matrices', ...
                         getfactorizationname(k));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            % with R columns
            if any(cellfun('size', U, 2) ~= size(U{1}, 2))
                adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                         ['Factorization %s: in a cpd all factors should have the ' ...
                          'same number (R) of columns (now: %s)'], ...
                         getfactorizationname(k), printsize(cellfun('size',U,2))); 
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            % and of appropriate size
            if ~sizeequals(getsize(U), getsize(model.factorizations{k}.data))
                adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                         ['Factorization %s: the size of the tensor generated by ' ...
                          'the factors %s should match the size of the data %s'], ...
                         getfactorizationname(k), printsize(getsize(U)), ...
                         printsize(getsize(model.factorizations{k}.data)));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
          case 'll1'
            L = model.factorizations{k}.L;
            U = cache.factors.expanded([model.factorizations{k}.factors{:}]);
            % factor matrices should be matrices
            if any(~cellfun(@ismatrix, U)) || any(~cellfun(@ isnumeric, U))
                adderror('factorizations', location, 'sdf:consistency:wrongTypeFactors', ...
                         'Factorization %s: in a LL1 all factors should be numerical matrices', ...
                         getfactorizationname(k));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            % there should be three factors
            if length(U) ~= 3
                adderror('factorizations', location, 'sdf:consistency:wrongNumberFactors', ...
                         ['Factorization %s: in a LL1 there should be exactly ' ...
                          'three factors'],...
                         getfactorizationname(k)); 
                cache.factorizations.errors{k}{end+1} = errorlist{end};
                continue;
            end
            % with sum(L) or length(L) columns
            if any(cellfun('size', U(1:2), 2) ~= sum(L)) || size(U{3},2) ~= length(L)
                adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                         ['Factorization %s: in a LL1 the first two factors ' ...
                          'should have sum(L) columns, and the third factor ' ...
                          'should have length(L) columns'],...
                         getfactorizationname(k)); 
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            % and of appropriate size
            if ~sizeequals(getsize(U), getsize(model.factorizations{k}.data))
                adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                         ['Factorization %s: the size of the tensor generated by ' ...
                          'the factors %s should match the size of the data %s'], ...
                         getfactorizationname(k), printsize(getsize(U)), ...
                         printsize(getsize(model.factorizations{k}.data)));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
          case {'btd', 'lmlra'}
            size_data = getsize(model.factorizations{k}.data);
            U = cellfun(@(i) cache.factors.expanded([i{:}]), ...
                        model.factorizations{k}.factors, 'UniformOutput', ...
                        false);
            if strcmp(model.factorizations{k}.type, 'lmlra') && length(U) ~= 1
                adderror('factorizations', location, 'sdf:consistency:wrongTypeFactors', ...
                         ['Factorization %s: in a lmlra, there should be only ' ...
                         'one term'], ...
                         getfactorizationname(k));                
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            
            for r = 1:length(U)
                size_fm = cellfun('size', U{r}(1:end-1), 2);
                % all factor matrices should be matrices
                if any(~cellfun(@ismatrix, U{r}(1:end-1))) || ...
                        any(~cellfun(@ isnumeric, U{r}))
                    adderror('factorizations', location, 'sdf:consistency:wrongTypeFactors', ...
                             ['Factorization %s: in a %s, in term r=%d, all ' ...
                              'factor matrices should be numerical matrices'],...
                             getfactorizationname(k), model.factorizations{k}.type, r);                
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                end
                % factor matrices and core should have compatible dimensions
                if ~sizeequals(size_fm, size(U{r}{end}))
                    adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                             ['Factorization %s: in a %s, for term r=%d, the ' ...
                              'column dimensions of each factor matrix %s should ' ...
                              'match the dimensions of the core tensor %s.'],...
                             getfactorizationname(k), model.factorizations{k}.type, r, ...
                             printsize(size_fm), printsize(size(U{r}{end})));
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                end
                % all terms should have same dimensions as the data
                if ~sizeequals(getsize(U(r)), size_data)
                    adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                             ['Factorization %s: the size of the tensor generated by ' ...
                              'the factors for term r = %d, %s, should match the ' ...
                              'size of the data %s'], ...
                             getfactorizationname(k), r, printsize(getsize(U(r))), ...
                             printsize(getsize(model.factorizations{k}.data)));
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                end
            end
          case {'regL0', 'regL1', 'regL2'}
            U = cache.factors.expanded([model.factorizations{k}.factors{:}]);
            if isfield(model.factorizations{k}, 'data')
                % number of datasets already checked
                % dataset and factor should have compatible sizes
                if any(~cellfun(@(u,d) all(size(u)==size(d)), U, ...
                                model.factorizations{k}.data))
                    adderror('factorizations', location, 'sdf:consistency:dimensionMismatch', ...
                             ['Factorization %s: the dimensions of one of the ' ...
                              'factors and its corresponding dataset are not ' ...
                              'equal'], ...
                             getfactorizationname(k));
                    cache.factorizations.errors{k}{end+1} = errorlist{end};
                end
            end
            % all values should be numerical
            if any(~cellfun(@isnumeric, U))  
                adderror('factorizations', location, 'sdf:consistency:wrongTypeFactors',...
                         ['Factorization %s: only numeric types can be ' ...
                          'regularized'], getfactorizationname(k));
                cache.factorizations.errors{k}{end+1} = errorlist{end};
            end
            % Check extra fields
            if strcmpi(model.factorizations{k}.type, 'regL0')
                if isfield(model.factorizations{k}, 'sigma')
                    sigma = model.factorizations{k}.sigma;
                    if ~isnumeric(sigma) && ~isa(sigma, 'function_handle')
                        adderror('factorizations', {k, 'sigma'}, 'sdf:consistency:invalidSigma',...
                                 ['Factorization %s: sigma should be either ' ...
                                  'numerical or a function handle.'], ...
                                 getfactorizationname(k));
                    elseif isa(sigma, 'function_handle')
                        U = cache.factors.expanded([model.factorizations{k}.factors{:}]);
                        try 
                            if ~isnumeric(sigma(1,U))
                                adderror('factorizations', {k, 'sigma'}, ...
                                         'sdf:consistency:invalidSigma',...
                                         ['Factorization %s: if sigma is a ' ...
                                          'function handle, its results should ' ...
                                          'be a scalar.'],getfactorizationname(k));
                            end
                        catch 
                            adderror('factorizations', {k, 'sigma'}, ...
                                     'sdf:consistency:invalidSigma',...
                                     ['Factorization %s: if sigma is a function ' ...
                                      'handle, it should accept two input ' ...
                                      'arguments: the iteration k and a cell ' ...
                                      'of factors U.'],getfactorizationname(k));
                        end
                    end
                end
            end
        end
        
        % Step 3: check for unknown fields
        knownfields = {'type', 'factors', 'data', 'datatype', 'isstructured', ...
                      'L', 'sigma', 'issymmetric', 'relweight', 'weight'};
        fn = fieldnames(model.factorizations{k});
        for l = 1:length(fn)
            if ~any(cellfun(@(t) strcmp(fn{l},t), knownfields))
                warning('sdf:consistency:unusedField', ...
                        'The field %s in %s is ignored.', ...
                        fn{l}, createlink(getsyntaxline('factorizations', k)));
            end
        end
    end

    %% Check variable and factor usage
    % check used factors
    usedfactors = cellfun(@(f) f.factors, model.factorizations, 'UniformOutput', ...
                          false);
    usedfactors = serialize(usedfactors);
    cache.factors.inuse = false(1, length(model.factors));
    cache.factors.inuse(usedfactors) = true;
    % check used variables
    usedvars = cellfun(@(f) cellfun(@(sf) sf{1}, f, 'UniformOutput', false), ...
                       model.factors(cache.factors.inuse), 'UniformOutput', ...
                       false);
    usedvars = cellfun(@(f) f(:), usedvars, 'UniformOutput', false);
    usedvars = cellfun(@(f,c) f(~c), usedvars, model.isconstant(cache.factors.inuse), ...
                       'UniformOutput', false);
    usedvars = cell2mat(cat(1, usedvars{:}));
    cache.variables.inuse = false(1, length(model.variables));
    cache.variables.inuse(usedvars) = true;

    %% Check if optimization variables are used
    if all(~cache.variables.inuse) || all(cellfun(@isempty, ...
                                                  model.variables(usedvars)))
        warning('sdf:consistency:noOptimizationVariables', ...
                ['The provided model only depends on constant factors or empty ' ...
                 'variables.']);
    end
    notinuse = find(~cache.variables.inuse);
    if ~isempty(notinuse)
        if length(notinuse) == 1
            plural = '';
            verb = 'is';
            vars =  getvariablename(notinuse);
        else 
            plural = 's';
            verb = 'are';
            vars = arrayfun(@getvariablename, notinuse, 'UniformOutput', ...
                           false);
            vars = join(vars, ', ');
        end
        warning('sdf:consistency:unusedVariable', ...
                'The variable%s %s %s unused.', plural, vars, verb);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 3: Display model summary (optional)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if options.Print 
        
        %% Display variables
        % not displayed
        
        %% Display factors
        ofmt = '%-20s %-20s %-40s\n';
        fprintf(ofmt, 'Factor', 'Size', 'Comment');
        fprintf([repmat('-', 1, 80) '\n']);
        for k = 1:length(model.factors)
            link = escapestring(viewexpansion(k));
            link = sprintf('fprintf([''%s'']);', link);
            % print links in red
            link = strrep(link, '$2$', ''']); fprintf(2, [''');   
            link = strrep(link, '$1$', '\n'']); fprintf([''');
            comments = {};
            
            if all(model.isconstant{k}(:))
                comments{end+1} = 'Constant';
            elseif any(model.isconstant{k}(:))
                comments{end+1} = 'Partially constant';
            end 
            if ~cache.factors.inuse(k)
                comments{end+1} = 'Unused';
            end
            if cache.factors.origin(k) > 0
                comments{end+1} = sprintf('Trans. %d', ...
                                          cache.factors.origin(k));
            elseif cache.factors.origin(k) == -2
                comments{end+1} = 'Implicit';
            end
            suberr = cache.factors.suberrors{k};
            ofmt1 = '%-20s %-20s %s';
            fprintf(ofmt1, createlink(getfactorname(k, true), link, 20), ...
                    printsize(size(cache.factors.expanded{k})), ...
                    join(comments, ','));
            if ~(isempty(cache.factors.errors{k}) && ...
                    all(cellfun(@(c) all(cellfun(@isempty,c)), suberr(:))))

                if ~isempty(comments), fprintf(','); end
                comments = {};
                
                if ~isempty(cache.factors.errors{k})
                    for l = 1:length(cache.factors.errors{k})
                        e = cache.factors.errors{k}{l};
                        errorlink = escapestring([e.message '\n']);
                        errorlink = sprintf('fprintf(2, [''%s'']);', errorlink);
                        switch e.identifier
                          case 'sdf:consistency:concatError'
                            comments{end+1} = createlink('Cannot concatenate subfactors', errorlink);
                          otherwise
                            comments{end+1} = 'Unknown';
                        end
                    end
                end
                
                for J = 1:numel(cache.factors.suberrors{k})
                    if ~isempty(cache.factors.suberrors{k}{J})
                        errs = ~cellfun(@isempty, cache.factors.suberrors{k}{J});
                        errs = cache.factors.suberrors{k}{J}(errs(:));
                        for l = 1:length(errs)
                            e = errs{l};
                            errorlink = escapestring([e.message '\n']);
                            errorlink = sprintf('fprintf(2, [''%s'']);', errorlink);
                            switch e.identifier
                              case 'sdf:consistency:structError'
                                comments{end+1} = createlink(['A transformation ' ...
                                                    'failed'], errorlink);
                              case 'sdf:consistency:unknownTransformation'
                                comments{end+1} = createlink(['Unknown ' ...
                                                    'transformation'], errorlink);
                              otherwise 
                                comments{end+1} = 'Unknown';
                            end
                        end
                    end
                end
                ofmt2 = '%s\n';
                fprintf(2, ofmt2, join(comments,','));
            else 
                fprintf('\n');
            end
        end 
        fprintf('\n');
        
        %% Display factorizations
        ofmt = '%-20s %-10s %-20s %-20s\n';
        fprintf(ofmt, 'Factorization', 'Type', 'Factors', 'Comments');
        fprintf([repmat('-', 1, 80) '\n']);
        for k = 1:length(model.factorizations)
            factorization = model.factorizations{k};
            comments = cell(1, length(cache.factorizations.errors{k}));
            for l = 1:length(cache.factorizations.errors{k})
                e = cache.factorizations.errors{k}{l};
                lonk = strrep(e.message, '"', ''' char(34) ''');
                lonk = strrep([lonk '\n'], sprintf('\n'), '\n');
                lonk = sprintf('fprintf(2, [''%s'']);', lonk);
                switch e.identifier
                  case 'sdf:consistency:wrongTypeFactors'
                    comments{l} = createlink('Factor error', lonk);
                  case 'sdf:consistency:dimensionMismatch'
                    comments{l} = createlink('Dimension mismatch', lonk);
                  otherwise
                    comments{l} = createlink('Error', lonk);
                end
            end
            comments = join(comments, ', ');

            ofmt1 = '%-20s %-10s %-20s ';
            ofmt2 = '%-20s\n';
            fprintf(ofmt1, ...
                    getfactorizationname(k), ...
                    factorization.type, ...
                    printfactors(factorization.factors, 20));
            fprintf(2, ofmt2, comments);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Finally, fix some things
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if constant, remove transformations
    for k = 1:length(model.factors)
        for l = 1:length(model.factors{k})
            if model.isconstant{k}(l)
                model.factors{k}{l} = cache.factors.sequence{k}{l}(end);
            end
        end
    end
    
    % remove transform field
    tmp = struct2cell(model);
    fields = fieldnames(model);
    mask = ~cellfun(@(s) strcmp(s, 'transform'), fields);
    model = cell2struct(tmp(mask), fields(mask));
    
catch e
    if ~strcmp(e.identifier, 'sdf:dummy'), rethrow(e); 
    else 
        throwAsCaller(e.cause{1}); 
    end   
end

if options.Expand 
    varargout = {cache.factors.expanded};
elseif options.Internal
    varargout{1} = model;
    varargout{2} = errorlist;
else
    varargout{1} = isconsistent;
end

    function str = join(data, link, gobble)
    if nargin >= 3 && gobble 
        data = cellfun(@strtrim, data, 'UniformOutput', false);
    end
    if isempty(data), str = ''; return; end
    if ~iscell(data), data = num2cell(data); end
    if isnumeric(data{1}), str = num2str(data{1});
    else str = data{1}; end
    for joinidx = 2:length(data)
        if isnumeric(data{joinidx})
            str = strcat(str, link, num2str(data{joinidx}));
        else
            str = strcat(str, link, data{joinidx});
        end
    end
    end

    function facts = dereffacts(factorizationindex, facts, level)
    if nargin < 3, level = {}; end
    for fidx = 1:length(facts)
        factorfound = false;
        location = [{k, model.factorizations{k}.type}, level];
        if iscell(facts{fidx})
            % recursive definiton needed for e.g. btd and lmlra
            facts{fidx} = dereffacts(factorizationindex, facts{fidx}, [level fidx]);
            factorfound = true; % otherwise an error is thrown earlier
        elseif ischar(facts{fidx}) && ~isfield(model.names, 'factors')
            factorfound = false;
        elseif ischar(facts{fidx})  
            idx = findname(facts{fidx}, model.names.factors);
            if isempty(idx) 
                factorfound = false;
            else 
                facts{fidx} = idx;
                factorfound = true;
            end
        elseif isscalar(facts{fidx}) && isnumeric(facts{fidx})
            idx = facts{fidx};
            if ~(round(idx) == idx && idx > 0 && idx <= length(model.factors)) ...
                    || isempty(model.factors{idx}) || ...
                    isemptyfactor(model.factors{idx})
                factorfound = false;
            else 
                factorfound = true;
            end
        else 
            adderror('factorizations', location, 'sdf:syntax:unknownFactor', ...
                     ['Factorization %s contains non-integer or non-string ' ...
                      'factor references'], ...
                     getfactorizationname(factorizationindex));
        end
        if ~factorfound
            idx = trymatchvariable(facts{fidx});
            if isnan(idx)
                varname = facts{fidx};
                if isnumeric(varname), varname = num2str(varname); end
                adderror('factorizations', location, 'sdf:syntax:unknownFactor', ...
                         'Factorization %s depends on the unknown factor %s', ...
                         getfactorizationname(factorizationindex), varname); 
            elseif ~isfinite(idx)
                varname = facts{fidx};
                if isnumeric(varname), varname = num2str(varname); end
                adderror('factorizations', location, 'sdf:syntax:unknownFactor', ...
                         ['Factorization %s depends on the unknown factor or ' ...
                          'variable %s'], ...
                         getfactorizationname(factorizationindex), varname); 
            else % found matching variable 
                if isfield(model.names, 'factors')
                    if isfield(model.names, 'variables')
                        fname = model.names.variables{idx};
                    else 
                        fname = sprintf('IF%d', idx);
                    end
                    position = length(model.factors)+1;
                    model.names.factors{position} = fname;
                else 
                    if ischar(facts{fidx})
                        if isempty(model.factors) || isempty(model.factors{idx}) || ...
                                isemptyfactor(model.factors{idx}) || ...
                                (cache.factors.origin(idx) == -2 && ...
                                 model.factors{idx}{1}{1} == idx)
                            position = idx;
                        else 
                            position = length(model.factors)+1;
                        end
                    else 
                        position = idx;
                    end
                    if isempty(model.factors) && isfield(model.names,'variables')
                        model.names.factors{position} = model.names.variables{idx};
                    end
                end
                
                % Create new factor and update reference
                model.factors{position} = {{idx}};
                facts{fidx} = position;
                
                % adjust constant and origins
                model.isconstant{position} = false;
                cache.factors.origin(position) = -2;
            end
        end
    end
    end

    function factors = fmtfactors(factors, type)
    if ischar(factors)
        factors = {factors};
    elseif ~iscell(factors)
        factors = num2cell(factors);
    end 
    if any(strcmp(type, {'btd', 'lmlra'}))
        if all(cellfun(@(r)isscalar(r)&&isnumeric(r), factors) | cellfun(@ischar, factors)) % lmlra
            factors = {factors};
        elseif any(~cellfun(@iscell, factors)) % convert ranges
            for i = find(~cellfun(@iscell, factors))
                factors{i} = num2cell(factors{i});
            end
        end
    end
    end

    function varidx = trymatchvariable(varname)
    % Returns nan if implicit factors are not allowed, otherwise returns inf
    % if no matching variable if found. If a matching variable is found, the
    % variable index is returned.
    varidx = inf;
    if ~model.options.ImplicitFactors, varidx = nan; return; end
    if ischar(varname)
        if isfield(model.names, 'variables')
            idx = findname(varname, model.names.variables);
            if ~isempty(idx)
                varidx = idx;
            end
        end
    elseif isnumeric(varname)
        if (round(varname) == varname && varname > 0 && varname <= ...
             length(model.variables)) && ~isempty(model.variables{idx})
            varidx = varname;
        end        
    end
    end

    function loc =  getsyntaxline(type, idx, varargin)
    if isfield(model.names, type)
        fieldref = sprintf('.%s', model.names.(type){idx});
    elseif ~isempty(idx)
        if strcmp(type, 'factorizations') && ...
                strcmp(factorizationtype, 'structarray')
            fieldref = sprintf('(%d)', idx);
        elseif strcmp(type, 'factors') && ...
                factorstructure.isarray
            fieldref = sprintf('(%d)', idx);
        else 
            fieldref = sprintf('{%d}', idx);
        end
    else 
        fieldref = '';
    end
    loc = sprintf('%s.%s%s', modelname, type, fieldref);
    for ri = 1:length(varargin)
        if ischar(varargin{ri})
            loc = sprintf('%s.%s', loc, varargin{ri});                
        elseif length(varargin{ri}) > 1                
            loc = sprintf('%s{%s}', loc, mat2str(varargin{ri}));
        else 
            append = sprintf('%d,', varargin{ri});
            loc = sprintf('%s{%s}', loc, append(1:end-1));
        end
    end
    end

    function adderror(varargin)
    if isnumeric(varargin{2}) || iscell(varargin{2})
        type = varargin{1};
        if iscell(varargin{2}), 
            idx = varargin{2}{1};
            remidx = varargin{2}(2:end);
        else 
            idx = varargin{2}; 
            remidx = {};
        end
        varargin = varargin(3:end);
    else 
        type = [];
        idx = 0;
    end
    
    indent = '';
    line = '';
    if ~isempty(type)
        loc = getsyntaxline(type, idx, remidx{:});
        line = sprintf('Error in line %s:\n\n', createlink(loc));
        indent = '    ';
    end
    if options.PrintErrors && ~options.StopOnError
        fprintf(2, ['\n' line indent varargin{2} '\n\n'], varargin{3:end});
    end
    err = MException(varargin{1}, [line indent varargin{2}], ...
                                          varargin{3:end});
    errorlist{end+1} = err;
        
    isconsistent = false;
    
    if options.StopOnError
        dummyerr = MException('sdf:dummy', varargin{1});
        dummyerr = addCause(dummyerr, err);
        throw(dummyerr);
    end
    
    end

    function name = getvariablename(varargin)
    name = getname('variables', varargin{:});
    end

    function name = getfactorname(varargin)
    name = getname('factors', varargin{:});
    end

    function name = gettransformname(varargin)
    name = getname('transform', varargin{:});
    end

    function name = getfactorizationname(varargin)
    name = getname('factorizations', varargin{:});
    end

    function name = getname(nametype, idx, noid)
    if nargin <= 2, noid = false; end
    if isfield(model.names, nametype)
        if noid, 
            fmt = '%s';
            name = sprintf(fmt, model.names.(nametype){idx});
        else 
            fmt = '%s (id = %d)'; 
            name = sprintf(fmt, model.names.(nametype){idx}, idx);
        end
    else
        name = num2str(idx);
    end
    end

    function numeric = checknumeric(data)
    if iscell(data)
        numeric = all(cellfun(@checknumeric, data));
    else 
        numeric = isnumeric(data);
    end
    end

    function name = getfunname(f)
    name = func2str(f);
    % remove arguments in case of anonymous functions:
    name = regexp(name, '(?:@\([^)]*\))?([^()]+)(?:\([^)]*\))?', 'tokens');
    name = name{1}{1};
    end

    function sz = printsize(sz) 
    sz = mat2str(sz);
    sz = strrep(sz, '  ', ' ');
    sz = strrep(sz, ' ', 'x');
    end

    function link = createlink(text, cmd, len)
    if nargin == 1, cmd = text; end
    if nargin <= 2, len = length(text); end
    if displaylinks
        link = sprintf('<a href="matlab: %s">%s</a>', cmd, text);
    else 
        link = text;
    end
    link = [link repmat(' ', 1, max(len-length(text), 0))];
    end

    function link = createdisplink(text, data, len)
    if iscell(data)
        cmd = cell2str(data);
    else
        cmd = ten2str(data);
    end 
    if nargin < 3
        link = createlink(text, cmd);
    else 
        link = createlink(text, cmd, len);
    end
    end

    function str = cell2str(data)
    for d = 1:numel(data)
        if iscell(data{d}), tmp = cell2str(data{d}); 
        else tmp = ten2str(data{d}); end
        if d == 1, str = tmp;
        else str = [str ', ' tmp]; end
    end
    str = sprintf('reshape({%s}, %s)', str, mat2str(size(data)));
    end

    function str = ten2str(data)
    str = sprintf('reshape(%s, %s)', mat2str(data(:)), ...
                  mat2str(size(data)));
    end

    function fact = printfactors(nb, len)
    if nargin < 2, len = 0; end
    if any(cellfun(@iscell, nb))
        fact = cellfun(@(n)printfactors(n,len), nb, 'UniformOutput', false);
        lengths = cellfun(@(s) len - (length(s) - length(strtrim(s))), fact);
        fact = ['{' join(fact, '},{', true) '}'];
        fact = [fact, repmat(' ', 1, max(len-sum(lengths)-length(lengths)*3+1, 0))];
    else 
        names = cellfun(@(n) getfactorname(n, true), nb, 'UniformOutput', ...
                        false);
        lengths = cellfun(@length, names);
        links = cellfun(@(n) viewexpansion(n), nb, 'UniformOutput', false);
        links = cellfun(@escapestring, links, 'UniformOutput', false);
        links = cellfun(@(s) sprintf('fprintf([''%s'']);', s), links, ...
                       'UniformOutput', false);
        links = strrep(links, '$2$', ''']); fprintf(2, [''');
        links = strrep(links, '$1$', '\n'']); fprintf([''');
        names = cellfun(@(n,l) createlink(n,l), names, links, 'UniformOutput', false);
        fact = join(names, ',');
        fact = [fact repmat(' ', 1, max(len-sum(lengths)-length(lengths)+1, 0))];
    end
    end


    function isequal = sizeequals(sz1, sz2)
    % Test if datasets have equal size regardless regardless of squeezed dimensions
    % sz1 = getsize(d1);
    % sz2 = getsize(d2);
    if length(sz1) > length(sz2)
        sz2 = [sz2 ones(1, length(sz1)-length(sz2))];
    elseif length(sz1) < length(sz2)
        sz1 = [sz1 ones(1, length(sz2)-length(sz1))];
    end
    isequal = all(sz1 == sz2);
    end

    function expansion = viewexpansion(factor)
    expansion = '\n';
    
    % get factor 
    if ischar(factor)
        fidx = find(cellfun(@(s) strcmp(s,factor), model.names.variables));
    else 
        fidx = factor;
        factor = getfactorname(fidx);
    end
    sequence = cache.factors.sequence{fidx};
    viewlink = createdisplink('view', cache.factors.expanded{fidx});
    if numel(sequence) == 1
        expansion = sprintf('%sExpansion of factor %s (%s):\n', expansion, ...
                            factor, viewlink);
    else 
        expansion = sprintf('%sExpansion of factor %s with structure %s (%s):\n', expansion, factor, ...
                            printsize(size(sequence)), viewlink);
    end
    for si = 1:numel(sequence)
        subs = cell(1,max(2, ndims(sequence)));
        [subs{:}] = ind2sub(size(sequence), si);
        coord = join(repmat({'%d'},1,length(subs)),',');
        if numel(sequence{si}) == 1
            if model.isconstant{fidx}(si), 
                varname = '(Constant)'; 
            else 
                varname = getvariablename(model.factors{fidx}{si}{1}, true);
            end
            varname = createdisplink(varname, sequence{si}{1},15);
            expansion = sprintf(['%s  (' coord ')  %-15s %-10s\n'], expansion, subs{:}, varname, ...
                                printsize(size(sequence{si}{1})));
        else 
            for li = 1:numel(sequence{si})-1
                if li == 1
                    invarname = getvariablename(model.factors{fidx}{si}{1}, true);
                    prep = sprintf(['  (' coord ') '], subs{:});
                else 
                    invarname = 'temp';
                    prep = '        ';
                end
                if li < numel(sequence{si}) - 1;
                    outvarname = 'temp';
                else 
                    outvarname = 'Factor';
                end
                invarname = createdisplink(invarname, sequence{si}{li}, 15);
                outvarname = createdisplink(outvarname, sequence{si}{li+1}, 15);
                
                factcomments = '';
                if ~isempty(cache.factors.suberrors{fidx}{si}{li})
                    e = cache.factors.suberrors{fidx}{si}{li};
                    link = strrep(escapestring(['\n' e.message '\n']), '\', '\\');
                    link = sprintf('fprintf(2, [''%s'']);', link);
                    switch e.identifier
                      case 'sdf:consistency:structError'
                        factcomments = createlink('transformation failed', ...
                                                  link);
                      case 'sdf:consistency:unknownTransformation'
                        factcomments = createlink('Unknown transformation', link);
                      otherwise 
                        factcomments = 'Unknown error';
                    end
                    factcomments = ['$2$', factcomments, '$1$'];
                end
                
                expansion = sprintf('%s%s %-15s %-10s %-20s %-15s %-10s %s\n', ...
                                    expansion, prep, ...
                                    invarname, printsize(size(sequence{si}{li})), ... 
                                    getfunname(model.factors{fidx}{si}{li+1}), ...
                                    outvarname, printsize(size(sequence{si}{li+1})), ...
                                    factcomments);
            end
        end
    end
    if ~isempty(cache.factors.errors{fidx})
        switch cache.factors.errors{fidx}{1}.identifier
          case 'sdf:consistency:concatError'
            expansion = sprintf('%s$2$Failed to concatenate subfactors$1$\n', ...
                                expansion);
        end
    end
    end

end

function str = escapestring(str)
str = strrep(str, '''', '''''');
str = strrep(str, '"', ''' char(34) ''');
str = strrep(str, sprintf('\n'), '\n');
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

% Local Variables:  
% mode: matlab           
% matlab-indent-function-body: nil
% End:              
