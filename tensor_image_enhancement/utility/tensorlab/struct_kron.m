function [x,state] = struct_kron(z,task)
%STRUCT_KRON Kronecker-product of two or more matrices.
%   [x,state] = STRUCT_KRON(z,[]) computes x as the Kronecker product
%   kron(z{1},z{2}). If length(z) = N > 2, x is computed as
%   kron(kron(...(kron(z{1},z{2}),z{3}),...),z{N}). The structure state stores
%   information which is reused in computing the right and left Jacobian-vector
%   products. The variables need to be vectors.
%
%   STRUCT_KRON(z,task) computes the right or left Jacobian-vector product
%   of this transformation, depending on the structure task. Use the
%   structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Otto Debals         (Otto.Debals@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

% Version History:
% - 2016/03/15   NV      Extended to N matrices
% - 2014/03/25   OD      Initial version for two vectors

if nargin < 2, task = []; end
if length(z) < 2
    error('struct_kron:z', 'z should be a cell with at least two matrices');
end

state = [];
sz1 = cellfun('size', z, 1);
sz2 = cellfun('size', z, 2);

right = isstruct(task) && isfield(task,'r') && ~isempty(task.r);
left  = isstruct(task) && isfield(task,'l') && ~isempty(task.l);

if ~left && ~right
    x = kron(z{1},z{2});
    for n = 3:length(z)
        x = kron(x, z{n});
    end
elseif right
    x = zeros(prod(sz1),prod(sz2));
    for n = 1:length(z)
        tmp = z;
        tmp{n} = task.r{n};
        xtmp = kron(tmp{1},tmp{2});
        for m = 3:length(z)
            xtmp = kron(xtmp,tmp{m});
        end
        x = x + xtmp;
    end
elseif left
    permdim = length(z):-1:1;
    permdim = [permdim; permdim+length(z)];
    task.l = reshape(task.l, [sz1(end:-1:1) sz2(end:-1:1)]);
    task.l = permute(task.l, permdim(:));
    task.l = reshape(task.l, sz1.*sz2);
    x = cell(1, length(z));
    z = cellfun(@(z) z(:), z, 'UniformOutput', false);
    for n = 1:length(z)
        x{n} = mtkronprod(task.l, z, n);
        x{n} = reshape(x{n}, sz1(n), sz2(n));
    end
end

end