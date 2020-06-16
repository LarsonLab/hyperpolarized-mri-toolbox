function [x,state] = struct_kr(z,task)
%STRUCT_KR Khatri-Rao-product of two or more matrices.
%   [x,state] = STRUCT_KR(z,[]) computes x as the Khatri-Rao product kr(z{:}) of
%   N = length(z) matrices with the same number of columns R. The structure
%   state stores information which is reused in computing the right and left
%   Jacobian-vector products. The variables need to be vectors.
%
%   STRUCT_KR(z,task) computes the right or left Jacobian-vector product
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

%   Authors:  Martijn Boussé      (Martijn.Bousse@esat.kuleuven.be)
%             Frederik Van Eeghem (Frederik.VanEeghem@kuleuven.be)
%             Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%             Otto Debals         (Otto.Debals@esat.kuleuven.be)
%             Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

% Version History:
% - 2016/03/15   NV      Simplification and speed up code + error tests
% - 2015/02/23   FVE     Extension to N > 2 complex matrices
% - 2015/02/18   MB      Initial version

if nargin < 2, task = []; end

if length(z) < 2
    error('struct_kr:z','z should be a cell with at least two matrices');
end
R = size(z{1},2);
if any(cellfun('size', z, 2)~=R)
    error('struct_kr:z','All matrices in the cell z should have R columns');
end

state = [];
right = isstruct(task) && isfield(task,'r') && ~isempty(task.r);
left  = isstruct(task) && isfield(task,'l') && ~isempty(task.l);

if ~left && ~right
    x = kr(z{:});
elseif right
    x = kr(task.r{1},z{2:end});
    for n = 2:length(z)
        tmp = z;
        tmp{n} = task.r{n};
        x = x + kr(tmp{:});
    end
elseif left
    sz = cellfun('size', z, 1);
    task.l = reshape(task.l, [sz(end:-1:1),R]);
    x = cell(1,length(z));
    z = [z(end:-1:1), eye(R)];
    for n = 1:length(z)-1
        x{n} = mtkrprod(task.l, z, length(z)-n);
        x{n} = reshape(x{n}, sz(n), R);
    end
end

end