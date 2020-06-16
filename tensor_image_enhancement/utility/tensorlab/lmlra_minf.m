function [U,S,output] = lmlra_minf(T,U0,S0,varargin)
%LMLRA_MINF LMLRA by unconstrained nonlinear optimization.
%   [U,S,output] = lmlra_minf(T,U0,S0) computes the factor matrices U{1},
%   ..., U{N} and core tensor S belonging to a low multilinear rank
%   approximation of the N-th order tensor T by minimizing 
%   0.5*frob(T-lmlragen(U,S))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the factor matrices U0{n} and core tensor S0. The
%   structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   lmlra_minf(T,U0,S0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%      options.Normalize =     - Normalize the result (factor matrices
%        [{true}|false]          orthogonal)
% 
%   See also lmlra_nls.

%   Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.
%
% Version History:
% - 2016/02/10   NV      Extracted core

    
    
% Check the options structure
p = inputParser;
p.addOptional('Algorithm', @minf_lbfgsdl);
p.addOptional('M', nan);
p.addOptional('MaxIter', 500);
p.addOptional('TolFun', 1e-12);
p.addOptional('TolX', 1e-8);
p.parse(varargin{:});

fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);
options.OptimizationType = 'minf';
    
options = [fieldnames(options)'; struct2cell(options)'];
[U,S,output] = lmlra_core(T,U0,S0,options{:});

end


