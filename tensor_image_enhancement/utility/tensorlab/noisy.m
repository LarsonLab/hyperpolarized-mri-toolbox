function [Y,N] = noisy(X,SNR,dist,combined,blocksize,variance)
%NOISY Generate a noisy version of a given array.
%   [Y,N] = NOISY(X,SNR) computes a noisy version of X as Y = X + N, where
%   the noise term N is generated as sigma*randn(size(X)) if X is real and
%   as sigma*(randn(size(X))+randn(size(X))*1i) if X is complex. The scalar
%   sigma is chosen such that 10*log10((X(:)'*X(:))/(N(:)'*N(:))) = SNR dB.
%   By default, SNR is 20. If X is a cell array, a noisy version of each of
%   its elements is recursively computed and returned in the cell array Y.
%
%   NOISY(X,SNR,DIST) uses the distribution DIST as noise source. The noisy
%   is given by sigma*dist(size(X)) if X is real and sigma*(dist(size(X)) +
%   1i*dist(size(X))) if X is complex. DIST is a function handle accepting
%   a size vector as argument. Default is @randn.
%
%   NOISY(X,SNR,DIST,COMBINED) computes the real and imaginary parts of the
%   noise in one function call if combined is true, i.e., if true, the
%   noise is given by sigma*dist(size(X)). If false, the noise is given by
%   sigma*dist(size(X)) if X is real and by
%   sigma*(dist(size(X))+1i*dist(size(X))) if X is complex . Default is
%   false.
%
%   X = NOISY(X,SNR,DIST,COMBINED,BLOCKSIZE,VARIANCE) uses less memory to
%   add the noise. Instead of computing all of the noise elements at once
%   (by default), noise is added in blocks of length BLOCKSIZE, e.g., 1e5.
%   To scale the noise, i.i.d. noise with zero mean is assumed, and the
%   variance of the noise is given by VARIANCE. For Gaussian distributed
%   noise (randn) the variance is 1 (default), for other distributions the
%   variance has to be provided. The output argument N is not provided when
%   this implementation is chosen.
%
%   See also randn, rand.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/08   NV      Large scale variant and distribution added
    

% Check the options structure.
if nargin < 2, SNR = 20; end
if nargin < 3 || isempty(dist), dist = @randn; end
if nargin < 4 || isempty(combined), combined = false; end
if nargin < 5, blocksize = nan; end
if nargin < 6, variance = 1; end

if nargout > 1 && nargin >= 5
    error('noisy:largescale', ['The noise N is not an output argument if ' ...
                        'blocksize is given']);
end

% Add noise recursively.
if nargout == 1
    Y = addnoise(X);
else 
    [Y,N] = addnoise(X);
end

function [X,N] = addnoise(X)
    if iscell(X)
        N = X;
        for i = 1:length(X)
            [X{i},N{i}] = addnoise(X{i});
        end
        return; 
    end
    if isstruct(X), size_x = size(X.val); real = isreal(X.val);
    else size_x = size(X); real = isreal(X); end
    if isnan(blocksize)
        tmp = dist(size_x);
        if ~combined && ~real, tmp = tmp + dist(size_x)*1i; end
        if isstruct(X), N = X; N.val = tmp;
        else N = tmp; end
        scale = frob(X)*10^(-SNR/20)/frob(N);
        if isstruct(X), 
            N.val = N.val*scale;
            X.val = X.val+N.val;
        else 
            N = N*scale;
            X = X+N; 
        end
    else
        if isstruct(X), numelX = length(X.val); 
        else numelX = numel(X); end
        % The expected value of zero mean independent noise is sqrt(sum(variance))
        scale = frob(X)/sqrt(numelX*variance)*10^(-SNR/20);
        if ~combined && ~real, scale = scale / sqrt(2); end
        for k = blocksize:blocksize:numelX
            tmp = scale * dist([1,blocksize]);
            if ~combined && ~real, tmp = tmp + 1i*scale*dist([1,blocksize]); end
            if isstruct(X), 
                X.val(k-blocksize+1:k) = X.val(k-blocksize+1:k) + tmp(:);
            else 
                X(k-blocksize+1:k) = X(k-blocksize+1:k) + tmp;
            end
        end
        if isempty(k), k = 0; end
        tmp = scale * dist([1,numelX-k]);
        if ~combined && ~real, tmp = tmp + 1i*scale*dist([1,numelX-k]); end
        if isstruct(X), 
            X.val(k+1:end) = X.val(k+1:end) + tmp(:);
        else 
            X(k+1:end) = X(k+1:end) + tmp;
        end
    end
end
end
