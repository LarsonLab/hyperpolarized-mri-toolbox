function [U,output] = btd_core(T,U0,varargin)
%BTD_CORE Computational core for block term decomposition.
%   BTD_CORE should not be called directly. Use BTD_MINF or BTD_NLS instead. 
%
%   See also btd_minf, btd_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Format the tensor T.
unstructuredtypes = {'full', 'incomplete', 'sparse'};
type = getstructure(T);
isstructured = ~any(strcmpi(type, unstructuredtypes));
if ~isstructured, 
    T = fmt(T,true); 
    type = getstructure(T);
end
isincomplete = strcmpi(type, 'incomplete');
issparse = strcmpi(type, 'sparse');

% Check the initial factors U0
U0 = U0(:).';
if any(~cellfun(@iscell, U0))
    error('btd_core:U0', 'U0 should be a cell of cells');
end
U0 = cellfun(@(u) u(:).', U0, 'UniformOutput', false);
if ~isvalidtensor(U0, false, 'btd')
    error('btd_core:U0', 'U0 is not a valid BTD initialization. See isvalidtensor(U0)');
end
size_tens = getsize(T);
sz = getsize(U0);
if length(sz) < length(size_tens)
    error('btd_core:U0', 'length(U0{1})-1 should be >= getorder(T)');
end
if any(sz(1:length(size_tens)) ~= size_tens & sz(1:length(size_tens))>0)
    error('btd_core:U0', 'size(U0{r}{n},1) should be size(T,n) for n=1:getorder(T).');
end
if any(sz(length(size_tens)+1:end)~=1) 
    error('btd_core:U0', 'size(U0{r}{n},1) should be 1 for n=getorder(T)+1:length(U0{r})-1');
end
size_tens = [size_tens, ones(1, length(sz)-length(size_tens))];
N = length(U0{1})-1;
R = length(U0);

% parse options
p = inputParser;
p.addOptional('OptimizationType', 'nls');
p.addOptional('Algorithm', @nls_gndl);
p.addOptional('CGMaxIter', 10);
p.addOptional('Display', 0);
p.addOptional('M', nan);
p.addOptional('TolLargeScale', 0.02);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-6);
p.addOptional('TolAbs', 0);
p.KeepUnmatched = true;
p.parse(varargin{:});
fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);
    
% Set absolute tolerances
try 
    cache.T2 = frob(T,'squared');
catch e
    if strcmpi(e.identifier, 'frob:notImplemented');
        error('btd_core:notImplemented', ...
              ['btd_core does not support the structured tensor type %s, yet. Use ' ...
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

% Check optimization algorithms
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
        warning('sdf:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
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
        warning('sdf:unknownAlgorithm', ['The %s is not known by BTD. Use ' ...
                            'at your own risk'], options.Algorithm);
    end
  otherwise
    error('sdf:unknownOptimizationType', ...
          ['The optimization type %s is unknown. Only ''nls'' and ''minf'' ' ...
           'are supported now'], options.OptimizationType);
end

% Select optimization subroutines.
usestate = false;
if strcmpi(options.OptimizationType, 'nls')
    usestate = true;
    dF.JHJx = @JHJx;
    dF.JHF = @grad;
    if isnan(options.M), options.M = 'block-Jacobi'; end
    switch options.M
      case 'block-Jacobi', dF.M = @M_blockJacobi;
      otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
    end
    state(U0,true);
else 
    dF = @grad;
    if isfield(options, 'M') && isnan(options.M)
        options = rmfield(options, 'M'); 
    end
end
    
% Call the optimization method.
cache.offset = cell2mat(cellfun(@(t)cellfun(@numel,t(:).'),U0(:).', ...
    'UniformOutput',false));
cache.offset = cumsum([1 cache.offset]);
[U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
output.Name = func2str(options.Algorithm);

if output.info == 4 && isstructured
    warning('btd_core:accuracy', ...
            ['Maximal numerical accuracy for structured tensors reached. The ' ...
             'result may be improved using ful(T) instead of T.']);
end

function state(z,firstrun)
    
    if nargin == 2 && firstrun
        % Store the fraction of known elements.
        if isincomplete
            cache.scale = length(T.val)./prod(T.size);
        end
    end
    
    % Cache the factor matrices' Gramians.
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    cache.UHU = ...
        arrayfun(@(i,n,j)z{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    
    % Optionally cache some results for the block-Jacobi preconditioner.
    if ischar(options.M) || isa(options.M,'function_handle')
        [idx,jdx] = ndgrid(1:R,1:N);
        UHU = cache.UHU;
        cache.invSKS = arrayfun( ...
            @(r,n)inv(mtkronprod(z{r}{end},UHU(r,:,r),n)* ...
            conj(reshape(permute(z{r}{end},[1:n-1 n+1:N n]), ...
            [],size(z{r}{end},n)))),idx,jdx,'UniformOutput',false);
        cache.invUHU = arrayfun( ...
            @(r,n)inv(UHU{r,n,r}),idx,jdx,'UniformOutput',false);
    end

end

function fval = objfun(z)
% BTD objective function.
    if isstructured
        fval = abs(0.5*cache.T2 - real(inprod(T, z)) + 0.5*frob(z,'squared'));
        return
    end

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
    
% BTD scaled conjugate cogradient.
    if usestate, state(z); end
    if ~isstructured, E = cache.residual; end
    offset = cache.offset;
    grad = nan(offset(end)-1,1);
    cnt = 1;
    for r = 1:length(z) 
        V = z{r}(1:N);
        S = conj(z{r}{end});
        for n = 1:N
            if isstructured
                tmp = full(mtkronprod(z,V,n)-mtkronprod(T,V,n))* ...
                      reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));                
            else
                tmp = full(mtkronprod(E,V,n))* ...
                      reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));
            end
            grad(offset(cnt):offset(cnt+1)-1) = tmp(:);
            cnt = cnt+1;
        end
        if isstructured
            tmp = full(mtkronprod(z,V,0)-mtkronprod(T,V,0));
        else
            tmp = full(mtkronprod(E,V,0));
        end
        grad(offset(cnt):offset(cnt+1)-1) = tmp;
        cnt = cnt+1;
    end
    
end

function y = JHJx(z,x)
    
    % BTD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    offset = cache.offset;
    UHU = cache.UHU;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    x = deserialize(x);
    XHU = arrayfun(@(i,n,j)x{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    y = nan(offset(end)-1,1);
    cnt = 1;
    for ri = 1:R
        
        % Factor matrices.
        for ni = 1:N
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
            cnt = cnt+1;
        end
        
        % Core tensor.
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
        cnt = cnt+1;
        
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
    invSKS = cache.invSKS;
    invUHU = cache.invUHU;
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
    
    % If incomplete, approximate the effect of missing entries.
    if isincomplete
        x = x/cache.scale;
    end
    
end

function x = deserialize(x)
    off = cache.offset; tmp = x; cnt = 1;
    x = cell(1,R);
    for r = 1:R
        x{r} = cell(1,N+1);
        for n = 1:N
            x{r}{n} = reshape(tmp(off(cnt):off(cnt+1)-1),size(U0{r}{n}));
            cnt = cnt+1;
        end
        x{r}{end} = reshape(tmp(off(cnt):off(cnt+1)-1),size(U0{r}{end}));
        cnt = cnt+1;
    end
end

end
