function T = lmlragen(U,S,varargin)
%LMLRAGEN Generate full tensor given a core tensor and factor matrices.
%   T = lmlragen(U,S) computes the tensor T as the mode-n tensor-matrix
%   product of the factor matrices U{n} with the core tensor S.
%
%   T = lmlragen(U,S,ind) computes only the elements of T corresponding to the
%   indices in ind. T has the same shape as the indices.
%
%   T = lmlragen(U,S,i,j,k,...) computes only the elements T(i,j,k,...) from the
%   tensor. The number of indices should match the order of the tensor given
%   by length(U). The colon operator can be used as a string, e.g., for a
%   third order tensor, lmlragen(U,S,i,':',k) can be used. 
%
%   See also btdgen, cpdgen, ll1gen, ttgen, lmlrares.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   Version history:
%   - 2016/01/13  NV   Added incomplete lmlragen

T = btdgen({[U(:).',S]},varargin{:});
