function X = kron(U,varargin)
%KRON Kronecker product.
%   kron(A,B) returns the Kronecker product of two matrices A and B, of
%   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
%   block matrix in which the (i,j)-th block is defined as A(i,j)*B.
%
%   kron(A,B,C,...) and kron({A B C ...}) compute a string of Kronecker 
%   products A x B x C x ..., where x denotes the Kronecker product.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Otto Debals         (Otto.Debals@cs.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version history:
% - 02/01/2016   OD      Variable number of input arguments      
% - 21/07/2015   OD      Efficient updates for vector inputs     

if nargin<2
    if isnumeric(U), X = U; return; end
    if iscell(U), X = kron(U{:}); return; end
elseif nargin>2
    X = kron(U,varargin{1});
    for i = 2:numel(varargin)
        X = kron(X,varargin{i});
    end
    return;
end

A = U;
B = varargin{1};

[I,J] = size(A);
[K,L] = size(B);

if ~issparse(A) && ~issparse(B)
    % Both matrices are dense.
    if J==1 && L==1
        X = B*A.';
        X = X(:);
    elseif I==1 && K==1
        X = B.'*A;
        X = X(:).';
    elseif J==1 && K==1
        X = A*B;
    elseif I==1 && L==1
        X = B*A;
    else
        A = reshape(A,[1 I 1 J]);
        B = reshape(B,[K 1 L 1]);
        X = reshape(bsxfun(@times,A,B),[I*K J*L]);
    end
else
    
    % One of the matrices is sparse.
    [ia,ja,sa] = find(A);
    [ib,jb,sb] = find(B);
    ix = bsxfun(@plus,K*(ia(:)-1).',ib(:));
    jx = bsxfun(@plus,L*(ja(:)-1).',jb(:));
    if islogical(sa) && islogical(sb)
        X = sparse(ix,jx,bsxfun(@and,sb(:),sa(:).'),I*K,J*L);
    else
        X = sparse(ix,jx,double(sb(:))*double(sa(:).'),I*K,J*L);
    end
    
end
