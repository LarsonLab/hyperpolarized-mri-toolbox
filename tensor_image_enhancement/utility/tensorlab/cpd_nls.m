function [U,output] = cpd_nls(T,U0,varargin)
%CPD_NLS CPD by nonlinear least squares.
%   [U,output] = CPD_NLS(T,U0) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T by minimizing 0.5*frob(T-cpdgen(U))^2. The algorithm is
%   initialized with the factor matrices U0{n}. The structure output
%   returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   CPD_NLS(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.LargeScale    - If true, the Gauss-Newton or Levenberg-
%      = sum(size(T))*R>1e2    Marquardt steps are computed using a
%                              preconditioned conjugate gradient algorithm.
%                              Otherwise, a direct solver is used.
%      options.M =           - The preconditioner to use when
%      [{'block-Jacobi'}|...   options.LargeScale is true.
%       'block-SSOR'|false]
%      options.PlaneSearch = - A function handle to the desired CPD plane
%      [{false},@cpd_eps]      search algorithm. Only applicable to
%                              algorithms with trust region globalization.
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX and
%                              options.PlaneSearchOptions. See also help
%                              [options.Algorithm].
%
%   See also cpd_minf.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Optimization-based
%       algorithms for tensor decompositions: canonical polyadic
%       decomposition, decomposition in rank-(Lr,Lr,1) terms and a new
%       generalization," SIAM J. Opt., 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Check the options structure.
    p = inputParser;
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('MaxIter', 200);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('JHasFullRank', false);
    p.addOptional('M', 'block-Jacobi');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    options.OptimizationType = 'nls';
    
    options = [fieldnames(options)'; struct2cell(options)'];
    [U, output] = cpd_core(T,U0,options{:});
end
