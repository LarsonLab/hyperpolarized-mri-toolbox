% TENSORLAB
% Version 3.0, 2016-03-28
%
% BLOCK TERM DECOMPOSITION
% Algorithms
%    btd_core     - Computational core for block term decomposition.
%    btd_minf     - BTD by unconstrained nonlinear optimization.
%    btd_nls      - BTD by nonlinear least squares.
% Initialization
%    btd_rnd      - Pseudorandom initialization for BTD.
% Utilities
%    btdgen       - Generate full tensor given a BTD.
%    btdres       - Residual of a BTD.
%    frobbtdres   - Frobenius norm of residual for a BTD.
%
% CANONICAL POLYADIC DECOMPOSITION
% Algorithms
%    cpd          - Canonical polyadic decomposition.
%    cpd_als      - CPD by alternating least squares.
%    cpd_core     - Computational routines for CPD decomposition.
%    cpd_minf     - CPD by unconstrained nonlinear optimization.
%    cpd_nls      - CPD by nonlinear least squares.
%    cpd_rbs      - CPD by randomized block sampling.
%    cpd3_sd      - CPD by simultaneous diagonalization.
%    cpd3_sgsd    - CPD by simultaneous generalized Schur decomposition.
% Initialization
%    cpd_gevd     - CPD by a generalized eigenvalue decomposition.
%    cpd_rnd      - Pseudorandom initialization for CPD.
% Line and plane search
%    cpd_aels     - CPD approximate enhanced line search.
%    cpd_els      - CPD exact line search.
%    cpd_eps      - CPD exact plane search.
%    cpd_lsb      - CPD line search by Bro.
% Utilities
%    cpd_crb      - Diagonal Cramér-Rao bound approximation for CPD.
%    cpderr       - Errors between factor matrices in a CPD.
%    cpdgen       - Generate full tensor given a polyadic decomposition.
%    cpdres       - Residual of a polyadic decomposition.
%    frobcpdres   - Frobenius norm of residual for polyadic decomposition.
%    rankest      - Estimate rank.
%
% COUPLED/SYMMETRIC CANONICAL POLYADIC DECOMPOSITION
% Algorithms
%    ccpd_core   - Computational routines for coupled/symmetric CPD
%                  decomposition.
%    ccpd_minf   - Coupled CPD by unconstrained nonlinear optimization.
%    ccpd_nls    - Coupled/symmetric CPD by nonlinear least squares.
%
% COMPLEX OPTIMIZATION
% Nonlinear least squares
%    nls_gncgs    - Nonlinear least squares by Gauss-Newton with CG-Steihaug.
%    nls_gndl     - Nonlinear least squares by Gauss-Newton with dogleg trust
%                   region.
%    nls_lm       - Nonlinear least squares by Levenberg-Marquardt.
%    nlsb_gndl    - Bound-constrained NLS by projected Gauss-Newton dogleg 
%                   trust region.
% Unconstrained nonlinear optimization
%    minf_lbfgs   - Minimize a function by L-BFGS with line search.
%    minf_lbfgsdl - Minimize a function by L-BFGS with dogleg trust region.
%    minf_ncg     - Minimize a function by nonlinear conjugate gradient.
%    minf_sr1cgs  - Minimize a function by SR1 with CG-Steihaug.
% Utilities
%    deriv        - Approximate gradient and Jacobian.
%    ls_mt        - Strong Wolfe line search by More-Thuente.
%    mpcg         - Modified preconditioned conjugate gradients method.
%
% LOW MULTILINEAR RANK APPROXIMATION
% Algorithms
%    lmlra        - Low multilinear rank approximation.
%    lmlra_core   - Computational core for low multilinear rank approximation. 
%    lmlra_hooi   - LMLRA by higher-order orthogonal iteration.
%    lmlra_minf   - LMLRA by unconstrained nonlinear optimization.
%    lmlra_nls    - LMLRA by nonlinear least squares.
%    lmlra3_dgn   - LMLRA by a differential-geometric Newton method.
%    lmlra3_rtr   - LMLRA by a Riemannian trust region method.
%    mlsvd        - (Truncated) multilinear singular value decomposition.
%    mlsvds       - Multilinear singular value decomposition for sparse tensors.
%    mlsvd_rsi    - Sequentially truncated MLSVD using randomized subspace
%                   iteration.
% Initialization
%    lmlra_aca    - LMLRA by adaptive cross-approximation.
%    lmlra_rnd    - Pseudorandom initialization for LMLRA.
% Utilities
%    lmlraerr     - Errors between factor matrices in a LMLRA.
%    lmlragen     - Generate full tensor given a core tensor and factor 
%                   matrices.
%    lmlrares     - Residual of a LMLRA.
%    froblmlrares - Frobenius norm of residual of a LMLRA.
%    mlrank       - Multilinear rank.
%    mlrankest    - Estimate multilinear rank.
%
% DECOMPOSITION IN MULTILINEAR RANK-(Lr,Lr,1) TERMS
% Algorithms
%    ll1          - Decomposition in LL1 terms.
%    ll1_core     - Computational routines for LL1 decomposition.
%    ll1_minf     - LL1 decomposition by nonlinear unconstrained
%                   optimization.
%    ll1_nls      - LL1 decomposition by nonlinear least squares.
% Initialization
%    ll1_gevd     - LL1 by generalized eigenvalue decomposition.
%    ll1_rnd      - Pseudorandom initialization for LL1 decomposition.
% Utilities
%    ll1convert   - Convert LL1 decomposition between CPD and BTD format.
%    ll1gen       - Generate full tensor given as a LL1 decomposition.
%    ll1res       - Residual for a LL1 decomposition.
%    frobll1res   - Frobenius norm of the residual of a LL1 decomposition.
%
% STRUCTURED DATA FUSION
% Language parser
%    sdf_check         - SDF language parser and syntax/consistency checker. 
% Algorithms
%    ccpd_core         - Computational routines for coupled/symmetric CPD
%                        decomposition.
%    ccpd_minf         - Coupled CPD by unconstrained nonlinear optimization.
%    ccpd_nls          - Coupled/symmetric CPD using nonlinear least squares.
%    sdf_core          - Computational core for structured data fusion.
%    sdf_minf          - Structured data fusion by unconstrained nonlinear 
%                        optimization.
%    sdf_nls           - Structured data fusion by nonlinear least squares.
% Structure
%    struct_abs        - Absolute value.
%    struct_band       - Band matrix.
%    struct_cell2mat   - Convert the contents of a cell array into a matrix.
%    struct_conj       - Complex conjugate.
%    struct_cauchy     - Cauchy matrix.
%    struct_const      - Keep parts of z constant.
%    struct_ctranspose - Complex conjugate transpose.
%    struct_diag       - Diagonal matrix.
%    struct_exp        - Matrix with columns as exponentials.
%    struct_fd         - Finite differences.
%    struct_gram       - Gramian matrix.
%    struct_hankel     - Hankel matrix.
%    struct_inv        - Matrix inverse.
%    struct_invsqrtm   - Matrix inverse square root.
%    struct_invtransp  - Matrix inverse transpose.
%    struct_kr         - Khatri-Rao-product of two or more matrices.
%    struct_kron       - Kronecker-product of two or more matrices.
%    struct_LL1        - Structure of third factor matrix in a LL1 decomposition.
%    struct_log        - Natural logarithm.
%    struct_matvec     - Matrix-vector and matrix-matrix product.
%    struct_nonneg     - Nonnegative array.
%    struct_nop        - No operation.
%    struct_normalize  - Normalize columns to unit norm.
%    struct_orth       - Rectangular matrix with orthonormal columns.
%    struct_plus       - Plus.
%    struct_poly       - Matrix with columns as polynomials.
%    struct_power      - Array power.
%    struct_prod       - Hadamard product.
%    struct_rational   - Matrix with columns as rational functions.
%    struct_rbf        - Matrix with columns as sums of Gaussian RBF kernels.
%    struct_select     - Select entry from cell variable z.
%    struct_sigmoid    - Constrain array elements to an interval.
%    struct_sqrt       - Square root.
%    struct_sum        - Sum of elements.
%    struct_times      - Times.
%    struct_toeplitz   - Toeplitz matrix.
%    struct_transpose  - Transpose.
%    struct_tridiag    - Tridiagonal matrix.
%    struct_tril       - Lower triangular matrix.
%    struct_triu       - Upper triangular matrix.
%    struct_vander     - Vandermonde matrix.
%
% TENSOR UTILITIES
% Structured tensors
%    detectstructure  - Detect structure in a tensor.
%    getstructure     - Determine the type of a tensor.
%    isvalidtensor    - Check if the representation of a tensor is correct.
% Products
%    contract         - Mode-n tensor vector contraction.
%    inprod           - Inner product of two tensors.
%    mtkronprod       - Matricized tensor Kronecker product.
%    mtkrprod         - Matricized tensor Khatri-Rao product.
%    outprod          - Outer vector/matrix/tensor product.
%    tmprod           - Mode-n tensor-matrix product.
% Utilities
%    fmt              - Format data set.
%    frob             - Frobenius norm.
%    ful              - Convert formatted data set to an array.
%    getorder         - Order of a tensor.
%    getsize          - Dimensions of a tensor.
%    noisy            - Generate a noisy version of a given array.
%    ttgen            - Generates full tensor from TT format.
%
% TENSORIZATION
% Deterministic tensorization
%    decimate         - Decimation of vectors, matrices or tensors.
%    hankelize        - Hankelization of vectors, matrices or tensors.
%    loewnerize       - Loewnerization of vectors, matrices or tensors.
%    mat2tens         - Tensorize a matrix.
%    segmentize       - Segmentation of vectors, matrices or tensors.
%    vec2tens         - Tensorize a vector.
% Deterministic detensorization
%    dedecimate       - Recover decimated signal(s).
%    dehankelize      - Recover the signal(s) from an (approximate) Hankel
%                       matrix/tensor.
%    deloewnerize     - Recover the signal(s) from an (approximate) Loewner
%                       matrix/tensor.
%    desegmentize     - Recover segmented signal(s).
%    tens2mat         - Matricize a tensor.
%    tens2vec         - Vectorize a tensor.
% Tensorization with statistics
%    dcov             - Covariance matrices along specific dimensions.
%    cum3             - Third-order cumulant tensor.
%    cum4             - Fourth-order cumulant tensor.
%    scov             - Shifted covariance matrices.
%    stcum4           - Fourth-order spatio-temporal cumulant tensor.
%    xcum4            - Fourth-order cross-cumulant tensor.
%
% UTILITIES
% Clustering
%    gap              - Optimal clustering based on the gap statistic.
%    kmeans           - Cluster multivariate data using the k-means++ algorithm.
% Polynomials          
%    genpolybasis     - Polynomial basis.
%    polymin          - Minimize a polynomial.
%    polymin2         - Minimize bivariate and real polyanalytic polynomials.
%    polyval2         - Evaluate bivariate and univariate polyanalytic 
%                       polynomials.
%    polysol2         - Solve a system of two bivariate polynomials.
%    ratmin           - Minimize a rational function.
%    ratmin2          - Minimize bivariate and real polyanalytic rational 
%                       functions.
%    transform_poly   - Transform a polynomial.
% Various
%    crandn           - Complex normally distributed pseudorandom numbers.
%    dotk             - Dot product in K-fold precision.
%    fixedanglevect   - Vectors with fixed angle.
%    gevd_bal         - Generalized eigenvalue decomposition with balancing.
%    kr               - Khatri-Rao product.
%    kron             - Kronecker product.
%    sumk             - Summation in K-fold precision.
% Visualization        
%    slice3           - Visualize a third-order tensor with slices.
%    spy3             - Visualize a third-order tensor's sparsity pattern.
%    surf3            - Visualize a third-order tensor with surfaces.
%    visualize        - Visualize a higher-order tensor.
%    voxel3           - Visualize a third-order tensor with voxels.

