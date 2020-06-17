function ncm = noise_correlation( im_noise )
%Calculates the noise correlation matrix from multi-channel noise-only 
%prescan. Expects only a single slice, single timepoint dataset
%from a multi-channel coil.
%
%INPUTS:
%im_noise: Noise-only prescan of the form [X Y num_coils]
%
%OUTPUTS:
%noise: Noise correlation matrix of size [X Y]. 

if ndims(im_noise) > 3
    warning('# of dimensions > 3, using 1st slice of Z dimension');
    im_noise = im_noise(:,:,1);
end

ncoils = size(im_noise,3);

tmp = reshape(im_noise,numel(im_noise)./ncoils,ncoils);
tmp = permute(tmp,[2 1]);
M = size(tmp,2);

ncm = (1/(M-1))*(tmp*tmp');

% figure,imagesc(abs(ncm));
% axis image off
% title('Noise Covariance Matrix');

end

