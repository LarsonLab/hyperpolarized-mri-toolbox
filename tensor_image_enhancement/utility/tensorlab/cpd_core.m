function [U,output] = cpd_core(T,U0,varargin)
%CPD_CORE Computational routines for CPD decomposition
%   CPD_CORE should not be called directly. Use CPD_MINF or CPD_NLS instead. 
%
%   See also cpd_minf, cpd_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Optimization-based
%       algorithms for tensor decompositions: canonical polyadic
%       decomposition, decomposition in rank-(Lr,Lr,1) terms and a new
%       generalization," SIAM J. Opt., 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

unstructuredtypes = {'full', 'incomplete', 'sparse'};
type = getstructure(T);
isstructured = ~any(strcmpi(type, unstructuredtypes));
if ~isstructured, 
    T = fmt(T,true); 
    type = getstructure(T);
end
size_tens = getsize(T);

isincomplete = strcmp('incomplete',type);
issparse = strcmp('sparse',type);
    
% Check the initial factor matrices U0.
U0 = U0(:).';
N = length(U0);
R = size(U0{1},2);
if any(cellfun('size',U0,2) ~= R)
    error('cpd_core:U0','size(U0{n},2) should be the same for all n.');
end
if length(size_tens) > N
    error('cpd_core:U0', 'length(U0) should be >= getorder(T)');
end
sz = cellfun('size', U0, 1);
if any(sz(1:length(size_tens)) ~= size_tens & sz(1:length(size_tens))>0)
    error('cpd_core:U0', 'size(U0{n},1) should be size(T,n) for all n');
end
if any(sz(length(size_tens)+1:end)~=1) 
    error('cpd_core:U0', 'size(U0{n},1) should be 1 for all n > getorder(T)');
end

% Check the options structure.
p = inputParser;
p.addOptional('OptimizationType', 'nls');
p.addOptional('M', nan);
p.addOptional('MaxIter', 200);
p.addOptional('CGMaxIter', 15);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-8);
p.addOptional('TolAbs', 0);
p.addOptional('LargeScale', sum(size_tens)*R>1e2);
p.addOptional('TolLargeScale', 0.02);
p.addOptional('JHasFullRank', false);
p.addOptional('Display', 0);
p.addOptional('UseCPDI', false);
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

% Line and plane search settings
isfunc = @(f)isa(f,'function_handle');
% Adapt line/plane search if it is a CPD line/plane search.
if isfield(options,'LineSearch') && isfunc(options.LineSearch) && ...
   ~isempty(strfind(func2str(options.LineSearch),'cpd_'))
    linesearch = options.LineSearch;
    options.LineSearch = @ls;
end
if isfield(options,'PlaneSearch') && isfunc(options.PlaneSearch) && ...
   ~isempty(strfind(func2str(options.PlaneSearch),'cpd_'))
    planesearch = options.PlaneSearch;
    options.PlaneSearch = @ps;
end
% Set absolute tolerances
try 
    cache.T2 = frob(T,'squared');
catch e
    if strcmpi(e.identifier, 'frob:notImplemented');
        error('cpd_core:notImplemented', ...
              ['cpd_core does not support the structured tensor type %s, yet. Use ' ...
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

% Process optimization type
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
        error('cpd_core:incompatibleAlgorithm', ...
              ['The %s method is incompatible with the nls optimization ' ...
               'type.'], func2str(options.Algorithm));
    elseif ~ismember(func2str(options.Algorithm), nlsfun)
        warning('cpd_core:unknownAlgorithm', ['The %s is not known by CPD. Use ' ...
                            'at your own risk'], options.Algorithm);
    end
  case 'minf'
    if ~isfield(options, 'Algorithm')
        options.Algorithm = @minf_lbfgsdl;
    elseif ismember(func2str(options.Algorithm), nlsfun)
        error('cpd_core:incompatibleAlgorithm', ...
              ['The %s method is incompatible with the minf optimization ' ...
               'type.'], func2str(options.Algorithm));
    elseif ~ismember(func2str(options.Algorithm), minffun)
        warning('cpd_core:unknownAlgorithm', ['The %s is not known by CPD. Use ' ...
                            'at your own risk'], options.Algorithm);
    end
  otherwise
    error('cpd_core:unknownOptimizationType', ...
          ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
           'are supported now'], options.OptimizationType);
end


% Call the optimization method.
usestate = false;
cache.offset = cumsum([1 cellfun(@numel,U0(:).')]);
if strcmpi(options.OptimizationType, 'nls')
    if options.UseCPDI
        if options.LargeScale, dF.JHJx = @JHJxi; 
        else dF.JHJ = @JHJi; end
    else
        if options.LargeScale, dF.JHJx = @JHJx; 
        else dF.JHJ = @JHJ; end
    end
   
    dF.JHF = @grad;
    if isnan(options.M), options.M = 'block-Jacobi'; end
    switch options.M
      case 'block-SSOR', dF.M = @M_blockSSOR;
      case 'block-Jacobi', dF.M = @M_blockJacobi;
      case 'Jacobi', dF.M = @M_Jacobi;
      otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
    end
    usestate = true;
    state(U0,true);
else 
    dF = @grad;
    if isfield(options, 'M') && isnan(options.M)
        options = rmfield(options, 'M'); 
    end
end 
[U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
output.Name = func2str(options.Algorithm);

if output.info == 4 && isstructured
    warning('cpd_core:accuracy', ...
            ['Maximal numerical accuracy for structured tensors reached. The ' ...
             'result may be improved using ful(T) instead of T.']);
end

% Computational core
function state(z,firstrun)

    if nargin == 2 && firstrun
        % Store the fraction of known elements.
        if isincomplete
            cache.scale = length(T.val)./prod(T.size);
        end
        if options.UseCPDI
            [cache.j, cache.i] = cellfun(@sort, T.sub, 'UniformOutput', 0);
            cache.i = cellfun(@double, cache.i, 'UniformOutput', false);
            cache.j = cellfun(@double, cache.j, 'UniformOutput', false);
        end
    end

    % Cache the factor matrices' Gramians.
    cache.UHU = zeros(N,R*R);
    for n = 1:N
        tmp = conj(z{n}'*z{n});
        cache.UHU(n,:) = tmp(:);
    end
    
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
    
    if options.UseCPDI
        J = cell(R,N);
        for n = 1:N
            dims = [1:n-1 n+1:N];
            i = cache.i{n};
            j = cache.j{n};
            for r = 1:R
                val = z{dims(1)}(T.sub{dims(1)},r);
                for k = dims(2:end)
                    val = val .* z{k}(T.sub{k},r);
                end
                val = val(i);
                tmp = sparse(i,j,val,length(T.val),size_tens(n));
                J{r,n} = tmp;                
            end
        end
        cache.J = cat(2, J{:});
    end 
    
end

function fval = objfun(z)
% CPD objective function.
    if isstructured
        fval = abs(0.5*cache.T2 - real(inprod(T, z)) + 0.5*frob(z,'squared'));
        return
    end
    
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = cpdgen(z); %z{1}*kr(z(end:-1:2)).';
        
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else
            fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:R
            tmp = z{1}(T.sub{1},r);
            for n = 2:length(z), tmp = tmp.*z{n}(T.sub{n},r); end
            fval = fval+tmp;
        end
    end
    if isincomplete
        E = T; E.val = fval;
        if ~isempty(E.matrix)
            E.matrix = sparse(double(E.sub{1}), ...
                double(1+idivide(E.ind-1,int64(size(E.matrix,1)))), ...
                double(fval),size(E.matrix,1),size(E.matrix,2));
        end
        cache.residual = E;
    else
        cache.residual = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function grad = grad(z)
    if usestate, state(z); end
    offset = cache.offset;    
    grad = nan(offset(end)-1,1);
    
    if isstructured
        for n = 1:length(z)
            tmp = mtkrprod(z,z,n) - mtkrprod(T,z,n);
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    else     
        % CPDn scaled conjugate cogradient.
        E = cache.residual;
        for n = 1:length(z)
            tmp = full(mtkrprod(E,z,n));
            grad(offset(n):offset(n+1)-1) = tmp(:);
        end
    end
end

function JHJ = JHJ(z)
    
    % CPD Jacobian's Gramian.
    UHU = conj(cache.UHU);
    JHJ = zeros(cache.offset(end)-1);
    for n = 1:N
        idxn = cache.offset(n):cache.offset(n+1)-1;
        Wn = reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]);
        JHJ(idxn,idxn) = kron(Wn,eye(size_tens(n)));
        for m = n+1:N
            idxm = cache.offset(m):cache.offset(m+1)-1;
            Wnm = reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),[R R]);
            JHJnm = bsxfun(@times,reshape(z{n},[size_tens(n) 1 1 R]), ...
                    reshape(conj(z{m}),[1 size_tens(m) R 1]));
            JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 1 R R]));
            JHJnm = permute(JHJnm,[1 3 2 4]);
            JHJnm = reshape(JHJnm,[size_tens(n)*R size_tens(m)*R]);
            JHJ(idxn,idxm) = JHJnm;
            JHJ(idxm,idxn) = JHJnm';
        end
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isincomplete
        JHJ = JHJ*cache.scale;
    end
    
end

function JHJ = JHJi(~)
    J = cache.J;
    JHJ = J'*J;
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
    
    % If incomplete, approximate the effect of missing entries.
    if isincomplete
        y = y*cache.scale;
    end
    
end

function y = JHJxi(~,x)
    J = cache.J;
    y = J*x;
    y = conj(y'*J);
    y = y(:);
end

function x = M_blockJacobi(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = nan(size(b));
    for n = 1:length(cache.offset)-1
        idx = cache.offset(n):cache.offset(n+1)-1;
        tmp = reshape(b(idx),[],size(cache.invW{1},1))*cache.invW{n};
        x(idx) = tmp(:);
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isincomplete
        x = x/cache.scale;
    end
    
end

function x = M_blockSSOR(z,b)
    
    % Solve Mx = b, where M is a block-Symmetric Successive Overrelaxation
    % preconditioner.
    % x = inv(D)*(U+D)*b
    B = cell(size(z));
    UHU = cache.UHU;
    BHU = nan(R,R,N);
    for n = 1:N
        B{n} = reshape(b(cache.offset(n):cache.offset(n+1)-1), size_tens(n), R);
        BHU(:,:,n) = B{n}'*z{n};
    end
    X = B;
    for n = 1:N-1
        Wsum = zeros(R);
        for m = n+1:N
            Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
            Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
        end
        Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
        X{n} = X{n}+(z{n}*Wsum)/Wn;
    end
    % x = (L+D)*x
    B = X;
    for n = 1:N
        Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
        BHU(:,:,n) = B{n}'*z{n};
        X{n} = B{n}*Wn;
    end
    for n = 2:N
        Wsum = zeros(R);
        for m = 1:n-1
            Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
            Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
        end
        X{n} = X{n}+z{n}*Wsum;
    end
    x = cell2mat(cellfun(@(x)x(:),X(:),'UniformOutput',false));
    
    % If incomplete, approximate the effect of missing entries.
    if isincomplete
        x = x/cache.scale;
    end
    
end

function [alpha,output] = ls(~,~,z,p,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = linesearch(T,z,p,state,options);
end

function [alpha,output] = ps(~,~,z,p,q,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = planesearch(T,z,p,q,state,options);
end

end
