function E = lmlrares(T,U,S,varargin)
%LMLRARES Residual of a LMLRA.
%   E = lmlrares(T,U,S) computes the residual tensor E as lmlragen(U,S)-T. If T
%   is an incomplete tensor, the residual is only computed for known elements.
%   If T is an efficient representation of a structured tensor, T is first
%   expanded using ful(T).
%
%   See also btdres, ll1res, cpdres, froblmlrares.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Compute the residual.
E = btdres(T,{[U(:).',S]},varargin{:});
