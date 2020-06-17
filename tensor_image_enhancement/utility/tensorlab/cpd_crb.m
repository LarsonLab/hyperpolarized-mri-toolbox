function [crb,varnoise] = cpd_crb(U, varnoise, varargin)
%CPD_CRB Diagonal Cramér-Rao bound approximation for CPD
%   C = CPD_CRB(U, varnoise) computes a diagonal approximation of the
%   Cramér-Rao bound under the assumption of Gaussian noise. U is a cell of N
%   factor matrices, each having R columns. The noise variance can be given
%   by varnoise.
%
%   [C, varnoise] = CPD_CRB(U, fval) computes the noise variance from the
%   function values fval as varnoise = mean(fval*2).
%
%   C = CPD_CRB(U, varnoise, V) matches the columns of the factors matrices
%   in U to the columns in the factor matrices in V by resolving the
%   permutation and scaling ambiguity (see cpderr). C is then computed using
%   the scaled and permuted factor matrices. 
%
%   CPD_CRB(U, varnoise, options) or CPD_CRB(U, varnoise, V, options) can be
%   used to set the following options:
%   - options.Method =              - use a (less accurate) diagonal
%        [{'diag'}, 'full']           approximation or the full gramian of
%                                     the Jacobian JHJ [1,3]. 
%   Key-value pairse can also be used to set the options.
% 
    
%   Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%              Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] N. Vervliet, L. De Lathauwer, "A randomized block sampling approach
%       to canonical polyadic decomposition of large-scale tensors," J. Sel.
%       Topics. Signal Process., IEEE, Vol. PP, No. 99, 2015.
%   [2] X. Liu and N. Sidiropoulos, "Cramér-Rao lower bounds for low-
%       rank decomposition of multidimensional arrays," IEEE Trans. Signal
%       Process., vol. 49, no. 9, pp. 2074-2086, Sept. 2001.
%   [3] P. Tichavský, A.-H. Phan, and Z. Koldovský, "Cramér-Rao-induced
%       bounds for CANDECOMP/PARAFAC tensor decomposition," IEEE
%       Trans. Signal Process., vol. 61, no. 8, pp. 1986-1997, Apr. 2013.
%
% Version History:
% - 2015/04/21   NV      Initial version
    
% parse options
    if nargin >= 3 && iscell(varargin{1}) && ...
                all(cellfun(@ismatrix, varargin{1}))
        % resolve permutation indeterminacy
        [~, P, D] = cpderr(varargin{1}, U);
        U = cellfun(@(u, d) u*P*d, U, D, 'UniformOutput', false);
        varargin = varargin(2:end);
    end

    p = inputParser;
    p.addOptional('Method', 'diag')
    p.parse(varargin{:});
    options = p.Results;
        
    % extract parameters from data 
    N = length(U);
    R = size(U{1},2);
    size_tens = cellfun(@(u) size(u,1), U);
    
    % determine noise level
    if ~isscalar(varnoise)
        varnoise = mean(varnoise*2);
    end

    offset = cumsum([1 cellfun(@numel, U)]);
    crb = zeros(offset(end)-1,1);
    if strcmpi(options.Method, 'diag')
        % We approximate JHJ by a diagonal, and invert only the diagonal
        nrms = cellfun(@(u) sum(abs(u).^2,1), U, 'UniformOutput', false);
        nrms = cat(1, nrms{:});
        for n = 1:N
            JHJd = prod(nrms([1:n-1 n+1:N],:),1);
            tmp = ones(size_tens(n),1)*(varnoise./JHJd);
            crb(offset(n):offset(n+1)-1) = tmp(:);
        end
    elseif strcmpi(options.Method, 'full')
        %% We compute the diagonal of the inverted JHJ matrix
        W = cellfun(@(u) u'*u, U, 'UniformOutput', false);
        W = cat(3, W{:});
        Ginv = cell(1,N);
        GinvZ = cell(1,N);
        Psi2 = cell(1,N);
        for n = 1:N
            tmp = prod(W(:,:,[1:n-1 n+1:N]),3);
            Ginv{n} = kron(inv(tmp), eye(size_tens(n)));
            Z = kron(eye(R), U{n});
            GinvZ{n} = kron(inv(tmp), U{n}); 
            Psi2{n} = Z'*GinvZ{n};
        end
        Psi = blkdiag(Psi2{:});
        P = sparse(reshape(1:R^2,R,R)',reshape(1:R^2,R,R), 1);
        K = zeros(N*R^2);
        for n = 1:N
            idxn = (n-1)*R^2 + (1:R^2);
            for m = n+1:N
                idxm = (m-1)*R^2 + (1:R^2);
                tmp = prod(W(:,:,[1:n-1 n+1:m-1 m+1:N]),3);
                K(idxn,idxm) = P*diag(tmp(:));
                K(idxm,idxn) = K(idxn,idxm);
            end
        end
        B = K*pinv(eye(N*R^2)+Psi*K);
        for n = 1:N
            idx = offset(n):offset(n+1)-1;
            idxb = (n-1)*R^2 + (1:R^2);
            tmp = GinvZ{n};
            crb(idx) = diag(Ginv{n}) - sum((tmp*B(idxb,idxb)).*tmp,2);
        end
        crb = crb*varnoise; 
    else 
        error('cpd_crb:unknownMethod', ['Currently, only the methods ''diag'' ' ...
                            'and ''full'' are known.'])
    end
end
