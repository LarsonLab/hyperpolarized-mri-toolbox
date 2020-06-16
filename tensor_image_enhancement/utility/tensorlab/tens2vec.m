function V = tens2vec(T,mode_row)
%TENS2VEC Vectorize a tensor
%   V = tens2vec(T,mode_row) vectorizes a full or sparse tensor T into a full or
%   sparse column vector V. The rows of V are obtained by looping over the
%   indices of T in the order mode_row. E.g., if A and B are two matrices and T
%   = cat(3,A,B), then tens2vec(T,1:3) is the vector [A(:);B(:)]. If
%   length(mode_row) < getorder(T), the remaining indices are added in the
%   original order. E.g., 
%       tens2vec(T,2) = tens2vec(T,[2 1 3]); % third-order tensor
%       tens2vec(T, [3 1]) = tens2vec(T,[3 1 2 4]); % fourth-order tensor
%
%   V = tens2vec(T) vectorizes a tensor T, where mode_row is chosen as the
%   sequence 1:getorder(T).
%
%   See also vec2tens, tens2mat.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/02   NV      Sparse tensors & less length(mode_row)<getorder(T)
% - 2014/02/01   LS      Initial version

if nargin < 2 && isnumeric(T)
    V = T(:);
elseif nargin < 2
    V = tens2mat(T,1:getorder(T));
    V = V(:);    
else 
    V = tens2mat(T,mode_row);
    V = V(:);
end
