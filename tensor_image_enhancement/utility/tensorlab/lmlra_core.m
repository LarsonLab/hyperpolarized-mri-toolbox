function [U,S,output] = lmlra_core(T,U0,S0,varargin)
%LMLRA_CORE Computational core for low multilinear rank approximation. 
%   LMLRA_CORE should not be called directly. Use LMLRA_MINF or LMLRA_NLS
%   instead. 
%
%   See also lmlra_minf, lmlra_nls.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.
%
% Version History:
% - 2016/02/10   NV      Extracted core and added U0 and S0 tests

% Parse options
p = inputParser;
p.addOptional('Normalize', true);
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;

% Normalize inputs
U0 = U0(:).';

% Check inputs
if any(~cellfun(@ismatrix, U0))
    error('lmlra_core:U0', 'U0{n} should be a matrix for all n');
end
if ~isnumeric(S0)
    error('lmlra_core:S0', 'S0 should be numeric');
end
size_tens = getsize(T);
sz = cellfun('size', U0, 1);
size_S = size(S0);
size_core = cellfun('size', U0, 2);
if length(sz) < length(size_tens)
    error('lmlra_core:U0', 'length(U0) should be >= getorder(T)');
end
if length(size_S) > length(size_core)
    error('lmlra_core:U0', 'length(U0) should be >= ndims(S0)');
end
if any(size_core(1:length(size_S)) ~= size_S)
    error('lmlra_core:U0', 'size(U0{n},2) should be size(S0,n) for all n');
end
if any(size_core(length(size_S)+1:end)~=1) 
    error('lmlra_core:U0', 'size(U0{n},2) should be 1 for all n > ndims(S0)');
end
if any(sz(1:length(size_tens)) ~= size_tens)
    error('lmlra_core:U0', 'size(U0{n},1) should be size(T,n) for all n');
end
if any(sz(length(size_tens)+1:end)~=1) 
    error('lmlra_core:U0', 'size(U0{n},1) should be 1 for all n > getorder(T)');
end

passopt = [fieldnames(p.Unmatched)'; struct2cell(p.Unmatched)'];
[U,output] = btd_core(T,{[U0,S0]}, passopt{:});
S = U{1}{end};
U = U{1}(1:end-1);

% normalize results
if options.Normalize
    [U,R] = cellfun(@(u) qr(u, 0), U, 'UniformOutput', false);
    S = lmlragen(R,S);
    [u,Sc] = mlsvd(S);
    
    sizeScflags = size(Sc)<size(S);
    if any(sizeScflags)
        ind = repmat({':'},1,ndims(S));
        for i = find(sizeScflags)
            ind{i} = 1:size(Sc,i);
            [u{i},~] = qr(u{i});
        end
        S = zeros(size(S)); S(ind{:}) = Sc;
    else S = Sc;
    end
    
    U(1:ndims(S)) = cellfun(@(u,v) u*v, U(1:ndims(S)), u, ...
                            'UniformOutput', false);
end 
end
