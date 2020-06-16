function [U,output] = ll1_rnd(size_tens, L, varargin)
%LL1_RND Pseudorandom initialization for LL1 decomposition.
%   U = ll1_rnd(size_tens, L) generates R = length(L) pseudorandom terms U{r} to
%   initialize algorithms that compute a (L,L,1)-decomposition of a third order
%   tensor. Each term U{r} is a cell array of three factor matrices U{r}{n} of
%   sizes size_tens(1) x L(r), size_tens(2) x L(r) and size_tens(3) x 1 followed
%   by a diagonal matrix U{r}{4} of size L(r).
%
%   ll1_rnd(T, L) is shorthand for ll1_rnd(getsize(T), L) if T is real. If T is
%   complex, then by default the terms will be generated using pseudorandom
%   complex numbers as well (cf. options).
%
%   ll1_rnd(size_tens, L, options) and ll1_rnd(T, L, options) may be used to set
%   the following options:
%
%      options.Real =         - The type of random number generator used to
%      [{@randn}|@rand|0]       generate the real part of the core tensor S
%                               and matrices U{n}. If 0, there is no real
%                               part.
%      options.Imag =         - The type of random number generator used to
%      [@randn|@rand|0|...      generate the imaginary part of the core
%       {'auto'}]               tensor S and matrices U{n}. If 0, there is
%                               no imaginary part. On 'auto', options.Imag
%                               is 0 unless the first argument is a complex
%                               tensor T, in which case it is equal to
%                               options.Real.
%      options.OutputFormat = - Format for the terms U: the BTD format returns
%      [{'btd'}|'cpd']          R=length(L) rank (L(r),L(r),1) terms with an
%                               identy core tensor, the CPD format returns 3
%                               factor matrices with sum(L), sum(L) and
%                               length(L) columns respectively.
%
%   See also cpd_rnd, lmlra_rnd, btd_rnd, btdgen.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

p = inputParser;
p.addOptional('OutputFormat', 'btd');
p.KeepUnmatched = true;
p.parse(varargin{:});
options = p.Results;
passopt = p.Unmatched;
    
isSizeVector = isnumeric(size_tens) && any(size(size_tens) == numel(size_tens));
if (isSizeVector && length(size_tens) ~= 3) || ...
        (~isSizeVector && getorder(size_tens) ~= 3)
    error('ll1_rnd:wrongOrder', 'll1_rnd only works for third order tensors');
end
size_core = arrayfun(@(l) [l, l, 1], L, 'UniformOutput', false);
[U, output] = btd_rnd(size_tens, size_core, passopt);
for r = 1:length(L), 
    U{r}{1} = U{r}{1}*U{r}{end};
    U{r}{end} = eye(L(r));
end
if strcmpi(options.OutputFormat, 'cpd')
    tmp = U;
    U = cell(1,3);
    for n = 1:3
        U{n} = cellfun(@(u) u{n}, tmp, 'UniformOutput', false);
        U{n} = cat(2, U{n}{:});
    end
elseif ~strcmpi(options.OutputFormat, 'btd')
    error('ll1_rnd:OutputFormat', 'Unknown output format %s', ...
          options.OutputFormat)
end
end
