function [U,output] = cpd_rnd(size_tens,R,varargin)
%CPD_RND Pseudorandom initialization for CPD.
%   U = cpd_rnd(size_tens,R) generates pseudorandom factor matrices U{1},
%   ..., U{N} of dimensions size_tens(n)-by-R that can be used to
%   initialize algorithms that compute the CPD of an N-th order tensor.
%
%   cpd_rnd(T,R) is shorthand for cpd_rnd(getsize(T),R) if T is real. If T is
%   complex, then U{n} will be generated as complex matrices by default (cf.
%   options).
%
%   cpd_rnd(size_tens,R,options) and cpd_rnd(T,R,options) may be used to
%   set the following options:
%
%      options.Real =        - The type of random number generator used to
%      [{@randn}|@rand|0]      generate the real part of each factor
%                              matrix. If 0, there is no real part.
%      options.Imag =        - The type of random number generator used to
%      [@randn|@rand|0|...     generate the imaginary part of each factor
%       {'auto'}]              matrix. If 0, there is no imaginary part.
%                              On 'auto', options.Imag is 0 unless the
%                              first argument is a complex tensor T, in
%                              which case it is equal to options.Real.
%      options.Orth =        - If true, the generated factor matrices are
%      [true|{false}|'auto']   orthogonalized using a QR factorization.
%                              The 'auto' option is deprecated. 
%      options.Angle = nan   - Fix the angle (in radians) between the factor
%                              matrices. This can be one angle for all factor
%                              matrices, or a vector with one angle for each
%                              factor matrix. A nan value will disable the
%                              fixed angle for a particular mode.
%      options.OptimalScaling
%        = [true, {false}]   - Optimally scale the factor vectors, such that they
%                              matches the input tensor as good as possible.
%                              This requires the first input to be a full,
%                              sparse or incomplete tensor.
%
%   See also btd_rnd, lmlra_rnd, cpdgen.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Version History:
% - 2016/03/18   NV      Default Orth = false
% - 2015/02/23   NV      Added optimal scaling and angles

% Process input.
isSizeVector = isnumeric(size_tens) && isvector(size_tens);
if ~isSizeVector
    T = size_tens;
    size_tens = getsize(T);
end
N = length(size_tens);

% Check inputs
isfunc = @(f)isa(f,'function_handle');
p = inputParser();
p.addOptional('Real', @randn);
p.addOptional('Imag', 'auto');
p.addOptional('Orth', 'auto');
p.addOptional('Angle', nan);
p.addOptional('OptimalScaling', false);
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;

% Check the options structure.
if ~isfunc(options.Real), options.Real = @zeros; end
if ischar(options.Imag) && strcmpi(options.Imag,'auto')
    if ~isSizeVector && ~isreal(ful(T,1))
        options.Imag = options.Real;
    else
        options.Imag = 0;
    end
end
if ischar(options.Orth) && strcmpi(options.Orth,'auto')
	options.Orth = false;
end
if options.OptimalScaling && isSizeVector, 
    warning('cpd_rnd:OptionDisabled', ...
            'Optimal scaling is disabled as no tensor is given.');
    options.OptimalScaling = false; 
end
if length(options.Angle) == 1, options.Angle = options.Angle*ones(1,N); end
if any(~isnan(options.Angle)), options.Orth = false; end

% Generate factor matrices.
if all(isnan(options.Angle))
    % Regular factor matrices
    U = arrayfun(@(n)options.Real(size_tens(n),R),1:N,'UniformOutput',0);
    if isfunc(options.Imag)
        Ui = arrayfun(@(n)options.Imag(size_tens(n),R),1:N,'UniformOutput',0);
        U = cellfun(@(ur,ui)ur+ui*1i,U,Ui,'UniformOutput',0);
    end
else 
    % Fixed angles factor matrices
    U = cell(1, N);
    for n = 1:N, 
        if ~isnan(options.Angle(n))
            U{n} = fixedAngleVect(size_tens(n), R, options.Angle(n), false); 
        else
            U{n} = options.Real(size_tens(n), R);
        end
    end
end 

% Scale factor matrices if requested
if options.OptimalScaling
    type = getstructure(T);
    if ~any(strcmpi(type, {'incomplete', 'sparse', 'full'}))
        error('cpd_rnd:OptimalScaling', ['Optimal scaling is currently only ' ...
                            'supports full, sparse and incomplete tensors.']);
    end
    if strcmpi(type, 'incomplete')
        A = U{1}(T.sub{1},:);
        for n = 2:length(U)
            A = A .* U{n}(T.sub{n},:);
        end
        b = T.val;
    elseif strcmpi(type, 'sparse') && 2*(sum(size_tens)*R)^2 < 1e9/8
        ind = randperm(prod(size_tens), 2*sum(size_tens)*R);
        sub = cell(1, length(size_tens));
        [sub{:}] = ind2sub(size_tens, ind);
        A = U{1}(sub{1},:);
        for n = 2:length(U)
            A = A .* U{n}(sub{n},:);
        end
        [~, ia, ib] = intersect(ind, T.ind);
        b = zeros(length(ind), 1);
        b(ia) = T.val(ib);
    elseif strcmpi(type, 'full') && numel(T)*sum(size_tens)*R < 1e9/8 % smaller than 1GB
        A = kr(U(end:-1:1));
        b = T(:);
    elseif strcmpi(type, 'full') && 2*(sum(size_tens)*R)^2 < 1e9/8
        ind = randperm(numel(T), 2*sum(size_tens)*R);
        sub = cell(1, length(size_tens));
        [sub{:}] = ind2sub(size_tens, ind);
        A = U{1}(sub{1},:);
        for n = 2:length(U)
            A = A .* U{n}(sub{n},:);
        end
        b = T(ind(:));
    else % Too large
        error('cpd_rnd:tooLarge', ['Tensor is too large to scale ' ...
                            'optimally']);
    end
    lambda = (A \ b).';
    U = cellfun(@(u) bsxfun(@times, sign(lambda).*nthroot(abs(lambda), N), u), U, ...
                'UniformOutput', false);
end

if any(~isnan(options.Angle))
    for n = 1:N
        rot = options.Real(size_tens(n));
        if isfunc(options.Imag), rot = rot + options.Imag(size_tens(n))*1i; end
        U{n} = orth(rot)*U{n};
    end
end

for n = 1:N*options.Orth
    if size(U{n},1) >= size(U{n},2), [U{n},~] = qr(U{n},0);
    else [Q,~] = qr(U{n}.',0); U{n} = Q.'; end
end
output = struct;
