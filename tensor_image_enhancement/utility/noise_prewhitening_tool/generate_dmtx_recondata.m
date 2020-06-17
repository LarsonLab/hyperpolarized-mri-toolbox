function [dmtx] = generate_dmtx_recondata(csi_phased)
% v1. generate noise correlation matrix for brain data
% use end of reconstructed spectra from edge voxels, last 3 timepoints
% ideally use raw data to get dmtx, but this should do for demo
size_csi_phased = size(csi_phased); % (x,1,y,f,chn,timepoint)
N_timepoints = 3; % timepoints to est noise correlation
Nt = size_csi_phased(6); % number of timepoints

% use edge voxels for noise
FID_noise = squeeze(csi_phased([1:2 end-1:end],1,[1:2 end-1:end],...
    :,:,Nt-N_timepoints+1:Nt));
FID_noise = permute(FID_noise,[3 1 2 5 4]);
FID_noise = reshape(FID_noise,[size_csi_phased(4) N_timepoints*16 size_csi_phased(5)]);

ncm = noise_correlation(FID_noise);
dmtx = ismrm_calculate_noise_decorrelation_mtx_from_covariance_mtx(ncm);

end