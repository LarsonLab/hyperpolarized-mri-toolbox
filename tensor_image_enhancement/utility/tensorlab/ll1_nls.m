function [U,output] = ll1_nls(T,U0,varargin)
%LL1_NLS LL1 decomposition by nonlinear least squares.
%   [U,output] = LL1_NLS(T,U0) computes R terms U{r} corresponding to a LL1
%   decomposition in BTD format of the third-order tensor T by minimizing
%   0.5*frob(T-ll1gen(U))^2. Each term U{r} is a cell array of N factor matrices
%   U{r}{n} with L(r) columns for n=1,2, and vector for n=3, followed by an
%   identity matrix U{r}{N+1}. The algorithm is initialized with the terms
%   matrices U0{r}. The structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   U = LL1_NLS(T, U0, L) computes the three factor matrices U{n}
%   corresponding to the LL1 decomposition in the CPD format of the
%   third-order tensor T. U{1} and U{2} both have sum(L) columns, U{3} has
%   length(L) columns. 
%
%   LL1_NLS(T,U0,options), LL1_NLS(T,U0,L,options), LL1_NLS(T,U0,'key', value)
%   and LL1_NLS(T,U0,L,'key',value)may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@nls_gndl}|...
%       @nls_gncgs|@nls_lm]
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX. See also
%                                help [options.Algorithm].
%      options.OutputFormat    - Either 'btd' or 'cpd'. If not given, the
%                                same format as the input is used. 
%
%   See also ll1_minf.

%   Authors: Nico Vervliet (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)


    if ~isempty(varargin) && isnumeric(varargin{1})
        [U, output] = ll1_core(T,U0,varargin{1}, 'Algorithm', @nls_gndl, ...
                               'OptimizationType', 'nls', varargin{2:end});
    else 
        [U, output] = ll1_core(T,U0,'Algorithm',@nls_gndl,'OptimizationType', ...
                               'nls', varargin{:});
    end
end