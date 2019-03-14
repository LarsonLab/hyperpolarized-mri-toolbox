function [matrix_X, k, diagv] = SVD_thres_kxy( matrix_X, mode, par,th)
% sigular value thresholding algrithm
% inputs: X, hankel matrix; k, k largest sigular values; lop to improve
% data consistence (need to double check)
% output: X, hankel matrix with rank of k
% initialized at Sep 9 2014 by Peng

% matrix_X = sparse(matrix_X);
% [U,S,V] = svds(matrix_X,par+1); %working version

%  [U,S,V] = rsvd(matrix_X,par+1); %not working
if  strcmp(mode, 'fixed_length') && (par == 1)
%     [U,S,V] = svds(matrix_X,par); 
    [U,S,V] = randomizedSVD(matrix_X,1, 2); %in TFOCS toolbox
else
    [U,S,V] = randomizedSVD(matrix_X,par+1, 2*par+2); %in TFOCS toolbox
end

%[U,S,V] = svdecon(matrix_X, par+1); %not work for sparse matrix
% [U,S,V] = svdsecon(matrix_X, par+1);

% [U,S,V] = svd(matrix_X, 'econ');

% opts.tol = 1e-8;
% opts.maxit = 150;
% [U,S,V] = lmsvd(matrix_X,par+1,opts);

clear matrix_X;

diagv = diag(S);
%plot(diagv(:));
if strcmp(mode, 'fixed_length')
    k = par;
elseif strcmp(mode, 'fixed_ratio')
    k=find(diagv<=(diagv(1)*par), 1);
elseif strcmp(mode, 'fixed_th')
%     th = 15;
    kfind=find(diagv<=th);
    k = min(kfind(:));
    if isempty(k)
        k = par;
    end
    if k < 1
        k = 1;
    end
    if k> par
        k = par;
    end
end
figure(222),
plot(diagv);

if strcmp(mode, 'fixed_th')
    diagv = diagv - th;
else
    if k > 1
        diagv = diagv - diagv(k);
    end
end

if k> 1
    diagv(k:end) = 0;
end

S = diag(diagv);



matrix_X = U*S*V';


end

