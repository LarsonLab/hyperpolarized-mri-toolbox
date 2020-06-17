function T = vec2tens(V,size_tens,mode_row)
%VEC2TENS Tensorize a vector.
%   T = vec2tens(V,size_tens,mode_row) tensorizes a full or sparse vector V into
%   a full or sparse tensor T of dimensions size_tens, given its vectorization
%   defined by mode_row. The rows of V should correspond to the elements of the
%   tensor by looping over the indices in the order mode_row (mode_col). E.g.,
%   if A and B are two matrices and V = [A(:);B(:)], then vec2tens(V,[size(A)
%   2],1:3) is the tensor T = cat(3,A,B).
%
%   T = vec2tens(V,size_tens) tensorizes a vector V, where mode_row is
%   chosen as the sequence 1:length(size_tens).
%
%   See also ten2vec, mat2tens.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/02   NV      Extension to sparse tensors
% - 2014/02/01   LS      Initial version

if nargin < 3 && ~issparse(V)
    T = reshape(V,size_tens);
else
    if nargin < 3
        mode_row = 1:length(size_tens);
    end 
    T = mat2tens(V,size_tens,mode_row);
end
