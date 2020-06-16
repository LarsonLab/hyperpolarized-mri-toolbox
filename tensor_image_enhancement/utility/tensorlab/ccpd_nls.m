function [U, output] = ccpd_nls(dataset, varargin)
%CCPD_NLS Coupled/symmetric CPD by nonlinear least squares.
%   U = CCPD_NLS(model) computes the coupled CPD defined by model. model is a
%   struct following the restricted SDF syntax: as factorization type only CPD
%   is allowed (but no CPDI), and no transformations, concatenations or constant
%   factors are allowed.
%   
%   U = CCPD_NLS(datasets, U0, cidx) computes the coupled CPD of the given
%   datasets (a cell of full, sparse, incomplete or structured tensors). Each
%   tensor is decomposed using a CPD with factors defined by cidx (a cell
%   with factor references with the same length as datasets), i.e., for
%   all n=1:length(datasets), datasets{n} is approximated by a CPD with
%   factors matricex U(i), i=cidx{n}. 
%
%   CCPD_NLS(model, options) and CCPD_NLS(datasets, U0, cidx, options) can be
%   used to set the following options (F is the number of datasets):
%    
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.PC = true     - Whether or not to use a preconditioner, if
%                              available, for computing the Gauss-Newton
%                              step.
%      options.RelWeights =  - By supplying relative weights, the weights
%      ones(1,F)               options.Weights are computed as follows:
%                              options.Weights(f) = options.RelWeights(f)/
%                              (sum(options.RelWeights)*numel(data_f)) for
%                              each of the F factorizations in the model.
%                              Relative weights are preferably set using the
%                              relweight field in a factorization.
%      options.Weights       - The weight of the F factorizations in the
%                              data fusion model. The SDF objective
%                              function is \sum_f 0.5*options.Weights(f)*
%                              frob(model_f-data_f)^2, where model_f and
%                              data_f are the fth factorization and
%                              corresponding tensor. By default, weights
%                              are provided by options.RelWeights, but
%                              options.Weights has precedence if supplied.
%                              Weights are preferably set using the weight
%                              field in a factorization.
%      options.IsSymmetric   - A logical array with length F, indicating that
%                              the dataset f has the same symmetry structure
%                              as indicated in the decomposition. (See
%                              example below.) IsSymmetric can also be set
%                              using the issymmetric field in a
%                              factorization. If nan, the algorithm checks if
%                              the dataset has the same symmetry as the
%                              decomposition. 
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX.
%                              See also help [options.Algorithm].
%
%   Example:
%
%   - Using SDF syntax
% 
%        model = struct;
%        model.variables = cpd_rnd([10 11 12 13], 5); 
%        model.factors = 1:4;
%        model.factorizations{1}.data = T1;
%        model.factorizations{1}.cpd  = 1:3;
%        model.factorizations{1}.issymmetric = true;
%        model.factorizations{2}.data = T2;
%        model.factorizations{2}.cpd  = [3 4 4];
%        model.factorizations{2}.issymmetric = false;
%
%        Ures = ccpd_nls(model);
%
%     The issymmetric field for factorization{1} indicates that the data
%     follows the same symmetry as the decomposition. As the decomposition
%     does not have symmetry, no symmetry is exploited. In factorization{2}
%     the decomposition is symmetric, but we explicitly state that the data
%     does not follow this symmetry, i.e., T(i,:,:) ~= T(i,:,:).' for some or
%     all i. 
%
%   - Using old syntax
%
%        data = {T1, T2};
%        U0   = cpd_rnd([10 11 12 13], 5);
%        cidx = {[1 2 3], [3 4 4]};
%       
%        Ures = ccpd_nls(data, U0, cidx, 'IsSymmetric', [true, false]);
 
% Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%          Otto Debals         (Otto.Debals@esat.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2016/01/10   NV      Initial version
    
    if nargin >= 2 && ~ischar(varargin{1}) && ~isstruct(varargin{1})
        if nargin < 3 || ischar(varargin{2}) || isstruct(varargin{2})
            error('ccpd_core:cidx', 'No coupling indices cidx provided');
        end 
        extras = varargin(1:2);
        varargin = varargin(3:end);
    else
        extras = {};
    end
    
    p = inputParser;
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('MaxIter', 200);
    p.addOptional('CGMaxIter', 15);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('TolLargeScale', 0.02);
    p.addOptional('JHasFullRank', false);
    p.addOptional('PC', true);
    p.addOptional('Display', 0);
    p.addOptional('M', 'block-Jacobi');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    options.OptimizationType = 'nls';
    
    options = [fieldnames(options)'; struct2cell(options)'];
    [U, output] = ccpd_core(dataset, extras{:}, options{:});

end