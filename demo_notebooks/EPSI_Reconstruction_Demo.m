% setup hyperpolarized-mri-toolbox

%IT WORKS, EXCCEPT AT LINE 98

%root_dir = sprintf('%s/',pwd);

% Setup data structures
fb_root_name = 'reconstruction/EPSI demo/csimage_in';
samp_pattern_name = 'hyperpolarized-mri-toolbox/reconstruction/EPSI demo/loc_samp_3d_dyn';
numReps = 18;

% setup additional parameters (from recon_cs3d_dyn_sivic_v1)

Navg = 1;

% compressed sensing weights for 
Itnlim = 16;
loopIter = 6;
TVWeight = 1e-6;
xfmWeight = 1e-2;
wavelet_time_scale = 2;  % similar between 0,1,2

plot_test = 0;

timedim = 4; % after squeeze

% load and setup data structures
load(samp_pattern_name); % File from blip generation

S = size(loc_samp_3d_dyn);
length_x = S(1); 
length_y = S(2); 
length_f = S(3);
length_z = 16;

% data_all = zeros(length_f, length_x, length_y, length_z, numReps);
mask_all = zeros(length_f, length_x, length_y, numReps);


for n = 1:numReps
    Isum = [1:Navg] + (n-1)*Navg;
    
    mask_all(:,:,:,n) = ...
        shiftdim(any(loc_samp_3d_dyn(:,:,:,Isum),4),2);
end

load(fb_root_name); % File from 3d undersample dataset

N = size(mask_all); % k-space and image size do not have to be the same
NN1 = 2^ceil(log(N(1))/log(2));
if numReps > 1
    NN4 = 2^ceil(log(N(4))/log(2)); % interpolate to diadic number for wavelets
else
    NN4 = 1;
end
NN = [NN1,N(2),N(3),NN4];

data_hybrid = fftshift(ifft(data_all, [], 4), 4);
kspace_hybrid = zeros(length_f, length_x, length_y, length_z, numReps);

% plot sampling pattern

repetition_plot = 8;
f_plot = 1:9;

clim = [0 1];
mask_plot = permute(mask_all(:,:,:,:),[2 3 4 1]); % X, Y, channel, F
%montage(mask_plot(:,:,repetition_plot,f_plot), 'DisplayRange',2)

%for f_plot = 1:5
%    subplot(1, 5, f_plot), imshow(permute(mask_all(f_plot,:,:,repetition_plot), [2 3 1]))
%end


spectra_plot_zerofilled = fftshift( fftshift( fftshift( fftshift(ifftn(ifftshift(ifftshift(ifftshift(data_all(:,:,:,:,repetition_plot), 2),3), 4)) ,1),2),3),4);

slice_plot = 9;
x_plot = 7;
y_plot = 9;

cplot(spectra_plot_zerofilled(:,x_plot,y_plot,slice_plot))
title('sample zero-filled complex spectrum')
figure()

total_carbon_image_zerofilled = squeeze(sum(abs(spectra_plot_zerofilled), 1));
imagesc(abs(total_carbon_image_zerofilled(:,:,9)))
%montage(permute(total_carbon_image_zerofilled, [1 2 4 3]))
title('sample zero-filled total carbon image')


%spectra_plot_zerofilled_slice = permute(abs(spectra_plot_zerofilled(:,:,:,slice_plot),[2 3 1]));

plot_voxels(permute(abs(spectra_plot_zerofilled(:,:,:,slice_plot)),[2 3 1]))

%%%%%%%%%%%%%
% L1 Recon
%%%%%%%%%%%%%

%kspace_all = recon_cs3d_dyn_sivic_v1(root_dir, fb_root_name, samp_pattern_name,numReps);

% % save dynamic data
% tmp_csreorder = read_one_image_ref('reconstruction/EPSI demo/dummy_ddf/dynth01_phased.cmplx');
% dynamic_cs.ddf = set_ddf_dyn_dim(tmp_csreorder.ddf, 'dynamic_cs', numReps);  
% dynamic_cs.ddf.specpoints = 59;
% dynamic_cs.img = kspace_all;
% dynamic_cs.img = fftshift(fft(dynamic_cs.img,[],1),1);
% for fftdim = 2:4
%     dynamic_cs.img = fftshift(ifft(dynamic_cs.img,[],fftdim),fftdim);
% end
% write_ddf_image_ex([root_dir 'dynamic_cs'],dynamic_cs);