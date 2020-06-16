function r = mlrank(T,tol)
%MLRANK Multilinear rank.
%   mlrank(A,tol) provides an estimate of the number of linearly independent
%   mode-n vectors in the full or sparse tensor T, for n = 1, ..., ndims(T).
%   More specifically, r(n) is the number of mode-n multilinear singular values
%   that are greater than tol(n)*max(sv{n}), where sv{n} are the mode-n
%   multilinear singular values of T. If tol is a scalar, it is set equal to
%   tol*ones(1,ndims(T)).
%
%   mlrank(A) uses the tolerance tol(n) = max(size(tens2mat(T,n)))*eps.
%
%   See also mlsvd, mlsvds, rank.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/02/07   NV    Added support sparse tensors

% Structured and incomplete tensors not supported yet.
type = getstructure(T);
if any(strcmpi(type, {'full', 'sparse', 'incomplete'}))
    T = fmt(T,true);
    if isstruct(T) && T.incomplete
        error('mlrank:T', ['Currently, this method only supports full and ' ...
                           'sparse tensors.']);
    end
else 
    error('mlrank:T', 'Currently, this method only supports full and sparse tensors.');
end

% Compute the multilinear singular values of T.
if isstruct(T)
    [~,~,sv] = mlsvds(T, getsize(T));
else
    [~,~,sv] = mlsvd(T);
end

% Set the tolerance.
sz = getsize(T);
N = getorder(T);
if nargin < 2, 
    tol = max(sz,prod(sz)./sz)*eps(class(ful(T,1))); 
end
tol = tol(:).';
if length(tol) == 1, tol = tol*ones(1,N); end
if length(tol) ~= N
    error('mlrank:tol', ...
          'tol should be a scalar, or a vector of length getorder(T).');
end

% Compute the multilinear rank.
r = arrayfun(@(n)sum(sv{n} > tol(n)*max(sv{n})),1:N);
