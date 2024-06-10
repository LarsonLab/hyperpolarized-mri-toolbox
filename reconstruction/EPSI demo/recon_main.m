clear all
close all
clc

root_dir = sprintf('%s/',pwd);
fb_root_name = 'csimage_in';
samp_pattern_name = 'loc_samp_3d_dyn';
numReps = 18;

kspace_all = recon_cs3d_dyn_sivic_v1(root_dir, fb_root_name, samp_pattern_name, ...
    numReps);

% save dynamic data
tmp_csreorder = read_ddf_image('dummy_ddf/dynth01_phased');
dynamic_cs.ddf = set_ddf_dyn_dim(tmp_csreorder.ddf, 'dynamic_cs', numReps);  
dynamic_cs.ddf.specpoints = 59;
dynamic_cs.img = kspace_all;
dynamic_cs.img = fftshift(fft(dynamic_cs.img,[],1),1);
for fftdim = 2:4,
    dynamic_cs.img = fftshift(ifft(dynamic_cs.img,[],fftdim),fftdim);
end
write_ddf_image_ex([root_dir 'dynamic_cs'],dynamic_cs);