function [U,output] = btd_minf(T,U0,varargin)
%BTD_MINF BTD by unconstrained nonlinear optimization.
%   [U,output] = btd_minf(T,U0) computes R terms U{r} corresponding to a
%   block term decomposition of the N-th order tensor T by minimizing
%   0.5*frob(T-btdgen(U))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the terms matrices U0{r}. The structure output returns
%   additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   btd_minf(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%
%   See also btd_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

    [U, output] = btd_core(T,U0,'Algorithm',@minf_lbfgsdl,'OptimizationType', ...
                           'minf', varargin{:});
end
