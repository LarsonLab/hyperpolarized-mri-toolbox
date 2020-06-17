function [U,output] = cpd_minf(T,U0,varargin)
%CPD_MINF CPD by unconstrained nonlinear optimization.
%   [U,output] = cpd_minf(T,U0) computes the factor matrices U{1}, ...,
%   U{N} belonging to a canonical polyadic decomposition of the N-th order
%   tensor T by minimizing 0.5*frob(T-cpdgen(U))^2. The algorithm is
%   initialized with the factor matrices U0{n}. The structure output
%   returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   cpd_minf(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.LineSearch =    - A function handle to the desired line
%      [{'auto'},@ls_mt|...      search algorithm. If the line search
%       @cpd_aels|@cpd_els|...   method has the prefix 'cpd_', it is
%       @cpd_lsb]                modified to be compatible with the
%                                optimization algorithm. Only applicable to
%                                algorithms with line search globalization.
%      options.PlaneSearch =   - A function handle to the desired CPD plane
%      [{false},@cpd_eps]        search algorithm. Only applicable to
%                                algorithms with trust region
%                                globalization.
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%
%   See also cpd_nls.

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
p.addOptional('Algorithm', @minf_lbfgsdl);
p.addOptional('M', nan);
p.addOptional('MaxIter', 500);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-8);
p.KeepUnmatched = true;
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);
options.OptimizationType = 'minf';

options = [fieldnames(options)'; struct2cell(options)'];
[U, output] = cpd_core(T,U0,options{:});
  
end
