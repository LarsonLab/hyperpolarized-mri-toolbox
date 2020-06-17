function T = ll1gen(U,varargin)
%LL1GEN Generate full tensor given as a LL1 decomposition
%   T = ll1gen(U) computes the third order tensor T given in the BTD format as
%   the sum of the R LL1-terms U{r} which are computed as tmprod(U{r}{N+1},
%   U{r}(1:N),1:N), where N is the order of the tensor T.
%
%   T = ll1gen(U,L) computes the third order tensor T given in the CPD format as
%   the sum of the R LL1-terms op(U{1}(:,i)*U{2}(:,i).',U{3}(:,r)) in which op()
%   is the outer product, and i=sum(L(1:r-1))+1:sum(L(1:r)) for r = 1:R. L hence
%   is a vector of length R = size(U{3},2) and sum(L) = size(U{1},2) =
%   size(U{2},2).
%
%   T = ll1gen(U,ind) for U in the BTD format and T = ll1gen(U,L,ind) for U
%   in the CPD format compute only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = ll1gen(U,i,j,k) for U in the BTD format and T = ll1gen(U,L,i,j,k) for U
%   in the CPD format compute only the elements T(i,j,k) the tensor. The number
%   of indices should match the order of the tensor given by length(U{1}) or
%   length(U) for the BTD format or the CPD format, respectively. The colon
%   operator can be used as a string, e.g., for a third order tensor and U in
%   the BTD format, ll1gen(U,i,':',k) can be used.
%
%   See also cpdgen, ll1gen, lmlragen, btdgen, ttgen, ll1res.
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/05   NV      Initial version
    
    if all(cellfun(@isnumeric, U)) 
        if nargin < 2
            error('ll1gen:L', ...
                  'If U is given in the CPD format, L should be given');
        end
        L = varargin{1};
        varargin = varargin(2:end);
        if any(cellfun('size', U(:).', 2) ~= [sum(L), sum(L), length(L)])
            error('ll1gen:U', ['cellfun(''size'', U, 2) should be [sum(L), ' ...
                               'sum(L), length(L)])']);
        end
        expvec = arrayfun(@(n) ones(1, L(n))*n, 1:length(L), 'UniformOutput', false); 
        expvec = cat(2, expvec{:});
        U{3} = U{3}(:,expvec);
        T = cpdgen(U,varargin{:});
    else 
        T = btdgen(U,varargin{:});
    end
end
