function [sol,output] = sdf_nls(model,varargin)
%SDF_NLS Structured data fusion by nonlinear least squares.
%   [sol,output] = sdf_nls(model) solves the data fusion problem described
%   by the structure model and returns the solution as the structure sol.
%   The model describes a data fusion problem with three fields:
%
%      model.variables
%
%         A structure or cell array of initializations for the variables in
%         the data fusion problem. Each field or cell in model.variables 
%         represents a variable. A variable may be an array (such as
%         scalars, vectors, matrices and tensors) or a (nested) cell array
%         of arrays.
%
%      model.factors
%
%         A structure or cell array of factors. Each field or cell in
%         model.factors represents a factor. A factor is described by a
%         cell array of subfactors. After each subfactor has been
%         generated, the factor is constructed as the cell2mat of its
%         subfactors. A subfactor is a cell array in which the first
%         element is a reference to a variable, and the following elements
%         represent a sequence of transformations of that variable. If
%         model.variables is a structure, then a reference to a variable is
%         a string corresponding to the field name of the desired variable.
%         If model.variables is a cell array, then a reference to a
%         variable is the index of the cell containing the desired
%         variable. Transformations are supplied by functions which
%         describe their (linearized) behaviour. For example, the
%         transformation @struct_inv computes the factor as the matrix
%         inverse of the selected variable. All transformations must have
%         the function signature struct_mytrans(z,task). Instead of
%         supplying a reference to a variable, it is also possible to
%         supply a constant factor in the form of an array.
%         
%      model.factorizations
%         
%         A structure of data sets to jointly factorize. Each field in
%         model.factorizations represents a factorization. A factorization
%         is described by a structure containing two fields. The first
%         field has the field name 'data' and contains the (dense, sparse
%         or incomplete) array which is to be factorized. If the array
%         contains many zeros or a NaN, it is internally converted to a
%         sparse or incomplete tensor with fmt. The second field's 
%         field name designates the type of factorization to compute of the
%         data set. Currently, two types are supported:
%
%         1. Canonical polyadic decomposition: a field 'cpd' should contain
%            a cell array of references to factors in model.factors. The
%            nth reference in the cell array corresponds to the nth factor
%            matrix of the CPD.
%
%         2. Block term decomposition: a field 'btd' should contain a cell
%            array of terms. Each term is itself a cell array containing
%            references to factors in model.factors. The nth reference in a
%            term corresponds to the nth factor matrix of that term. The
%            (N+1)th reference, where N is the number of dimensions of the
%            data set, in a term corresponds to the core tensor of that
%            term.
%
%         Additionally, two types of regularization are available:
%
%         3. L2 regularization: a field 'regL2' should contain a cell array
%            of references to factors in model.factors. This appends a term
%            to the objective function of the form 0.5*norm(F(:)-D(:),2)^2,
%            where F(:) and D(:) are the serialized factors in regL2 and
%            the data field respectively. The data field may be omitted, in
%            which case it is an all-zero vector.
%
%         4. L1 regularization: a field 'regL1' should contain a cell array
%            of references to factors in model.factors. This appends a term
%            to the objective function which is a smooth approximation of
%            0.5*norm(F(:)-D(:),1), where F(:) and D(:) are the serialized
%            factors in regL1 and the data field respectively. The data
%            field may be omitted, in which case it is an all-zero vector.
%
%   The algorithm proceeds to compute the joint factorization of the data
%   sets provided by models.factorizations by minimizing the sum of square
%   magnitude residuals of these factorizations. The output sol contains
%   the optimized variables in sol.variables and the corresponding factors
%   in sol.factors.
%
%   The structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   sdf_nls(model,options) may be used to set the following options:
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
%      options.Weights       - The weight of the F factorizations in the
%                              data fusion model. The SDF objective
%                              function is \sum_f 0.5*options.Weights(f)*
%                              frob(model_f-data_f)^2, where model_f and
%                              data_f are the fth factorization and
%                              corresponding tensor. By default, weights
%                              are provided by options.RelWeights, but
%                              options.Weights has precedence if supplied.
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX.
%                              See also help [options.Algorithm].
%
%   See also sdf_minf.

%   Authors: Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%            Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       J. Sel. Topics Signal Process., IEEE, Vol. 9, No. 4, pp. 586-600,
%       2015.

% Check the options structure.
    p = inputParser;
    p.addOptional('Algorithm', @nls_gndl);
    p.addOptional('MaxIter', 500);
    p.addOptional('CGMaxIter', 15);
    p.addOptional('TolFun', 1e-12);
    p.addOptional('TolX', 1e-8);
    p.addOptional('TolLargeScale', 0.02);
    p.addOptional('JHasFullRank', false);
    p.addOptional('PC', true);
    p.addOptional('M', 'block-Jacobi');
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    
    fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
    data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
    options = cell2struct(data, fn);
    options.OptimizationType = 'nls';
    
    options = [fieldnames(options)'; struct2cell(options)'];
    [sol, output] = sdf_core(model,options{:});
end
