function [csi_whitened] = noise_prewhitening(csi,dmtx)
% v1. noise pre-whitening for brain 2D EPSI data 20170831
% -----------------------
% Assumes data is csi(x,z,y,f,coil,t)

fprintf('noise pre-whitening ... \n');

fsize = size(csi);

if length(fsize) < 6
    Nt = 1;
else
    Nt = fsize(6);
end

Nd = 4;  % assumes x,y,z,f data

% to kspace
k_csi = csi;
for i_dim = 1:Nd
	k_csi = 1/sqrt(size(k_csi,i_dim))*fftshift(fft(ifftshift(k_csi,i_dim),[],i_dim),i_dim);
end

% pre-whitening kernal
k_csi_whitened = zeros(fsize);
% csi_whitened = ismrm_apply_noise_decorrelation_mtx(csi, dmtx);
for i_freq = 1:fsize(4)
for i_time = 1:Nt
    k_csi_whitened(:,:,:,i_freq,:,i_time) = ismrm_apply_noise_decorrelation_mtx(squeeze(k_csi(:,:,:,i_freq,:,i_time)), dmtx);
end
end

% to xspace
csi_whitened = k_csi_whitened;
for i_dim = 1:Nd
    csi_whitened = sqrt(size(csi_whitened,i_dim))*fftshift(ifft(ifftshift(csi_whitened,i_dim),[],i_dim),i_dim);
end

% original scale (not necessary)
csi_whitened = csi_whitened / max(abs(csi_whitened(:))) * max(abs(csi(:)));

end
