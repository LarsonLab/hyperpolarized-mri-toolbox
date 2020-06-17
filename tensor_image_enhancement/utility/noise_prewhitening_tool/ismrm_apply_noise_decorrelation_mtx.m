function [out] = ismrm_apply_noise_decorrelation_mtx(in, dmtx)
%
%  [out] = ismrm_apply_noise_decorrelation_mtx(in, dmtx)
%
%  Applies noise decorrlation matrix dmtx to data
%
%  INPUT:
%    - input   [nx, ny, ... , coils]   : Input data, last dimension is
%                                        coils
%    - dmtx    [coils,coils]           : Decorrelation matrix, e.g. as
%                                        produced by ismrm_calculate_noise_decorrelation_mtx
%
%  OUTPUT:
%    - out     [nx, ny, ..., coilx]    : Data with pre-whitened noise
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

orig_size = size(in);
ncoils = size(dmtx,1);
nelements = prod(orig_size)/ncoils;

if ~(ncoils == orig_size(end)),
    error('Number of coils in decorrelation matrix does not match the data');
end

in = permute(reshape(in,nelements,ncoils),[2 1]);
out = reshape(permute(dmtx*in,[2 1]),orig_size);

return 
