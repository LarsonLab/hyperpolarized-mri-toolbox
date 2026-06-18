% copying over brainweb
clear; close all;
addpath('../');
addpath('../brainweb_clone');

%% PARAMETERS
% tissue
mask = double(load('brainweb_fuzzy.mat').im_mask);

k_trans = [1, 0.2, 0.2;
           3, 0.4, 0.4];

% pk model params
tr = 4;
n_t = 30;
flips = repmat([20; 30; 30;], 1, n_t) .* (pi/180);
r1 = [1/30, 1/25, 1/25]; % pyr, lac, bic
k = [0, 0.05, 0.03; % lac, in order [vasc, gm, wm]
     0, 0.02, 0.01]; % bic
t_arrival = [0, 0, 0]; % vasc, gm, wm
t_bolus = 8;

% input function and mz0
n_tissues = 3;
input_function = zeros(n_tissues, n_t);
for i_tissue = 1:n_tissues
    input_function(i_tissue, :) = realistic_input_function(n_t, tr, t_arrival(i_tissue), t_bolus);
end
mz0 = [input_function(1), input_function(1)*.5, input_function(1)*.5;
       0, input_function(1)*.01, input_function(1)*.01;
       0, input_function(1)*.005, input_function(1)*.005];

% mri
coil_lim = [0.4, 1.2];
augment_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5]);

sample_size = [16 16 8; 16 16 8; 16 16 8];
snr = [150 40 20];
output_size = [32 32 8; 32 32 8; 32 32 8];

%% RUNNING THE MODEL -----------------------------------------------------------
% tissue
brain = tissue_structure("brain", mask, ["vasc", "gm", "wm"], [100,100,100]);
brain = brain.create_k_trans_map(k_trans);

% pk model
[met_images, ~, ~, ~, met_images_no_ktrans] = pk_model.run_pk_model(mz0, r1, k, flips, tr, brain, input_function=input_function, plot=true, met_names=["pyr", "lac", "bic"]);

% mri
cell_augment_params = namedargs2cell(augment_params); % unpack the augmentation parameters
met_images_aug = mri_system.augment(met_images, cell_augment_params{:});
met_images_coil_lim = mri_system.apply_coil_lim(met_images_aug, coil_lim, brain.Mask);
met_images_mres = mri_system.make_met_images_multires(met_images_coil_lim, sample_size);
met_images_mres_noise = mri_system.add_rician_noise(met_images_mres, snr);

met_images_upsampled = mri_system.upsample_to_output_size(met_images_mres_noise, output_size);




%% DISPLAY ---------------------------------------------------------------------

% tissue
brain.plot_alpha_composite_image(slice=50, order=[2,3,1]); % plot vasculature last

slices = 10:10:90;
figure(Name='kTRANS');
imagescn(brain.K_trans_map(:,:,slices), [0 max(brain.K_trans_map(:,:,slices), [], 'all')], [1 numel(slices)]);
colormap fire;

k_trans_dwnszd = imresize3(imresize3(brain.K_trans_map, [16 16 8]), [32 32 8]);
slices = 1:size(k_trans_dwnszd, 3);
figure(Name='kTRANS downsized');
imagescn(k_trans_dwnszd(:,:,slices), [0 max(k_trans_dwnszd(:,:,slices), [], 'all')], [1 numel(slices)]);
colormap fire;

%% met_images
% no ktrans
slices = round(size(met_images_no_ktrans, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (unprocessed)');
imagescn(squeeze(met_images_no_ktrans(:,:,slices,1,time_pts)), [0, max(met_images_no_ktrans(:,:,slices,1,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (unprocessed)');
imagescn(squeeze(met_images_no_ktrans(:,:,slices,2,time_pts)), [0, max(met_images_no_ktrans(:,:,slices,2,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (unprocessed)');
imagescn(squeeze(met_images_no_ktrans(:,:,slices,3,time_pts)), [0, max(met_images_no_ktrans(:,:,slices,3,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

% ktrans
slices = round(size(met_images, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (ktrans)');
imagescn(squeeze(met_images(:,:,slices,1,time_pts)), [0, max(met_images(:,:,slices,1,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (ktrans)');
imagescn(squeeze(met_images(:,:,slices,2,time_pts)), [0, max(met_images(:,:,slices,2,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (ktrans)');
imagescn(squeeze(met_images(:,:,slices,3,time_pts)), [0, max(met_images(:,:,slices,3,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

%% met images of various mri steps

% augmentations
slices = round(size(met_images_aug, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (aug)');
imagescn(squeeze(met_images_aug(:,:,slices,1,time_pts)), [0, max(met_images_aug(:,:,slices,1,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (aug)');
imagescn(squeeze(met_images_aug(:,:,slices,2,time_pts)), [0, max(met_images_aug(:,:,slices,2,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (aug)');
imagescn(squeeze(met_images_aug(:,:,slices,3,time_pts)), [0, max(met_images_aug(:,:,slices,3,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

% coil limits
slices = round(size(met_images_coil_lim, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (coil lim)');
imagescn(squeeze(met_images_coil_lim(:,:,slices,1,time_pts)), [0, max(met_images_coil_lim(:,:,slices,1,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (coil lim)');
imagescn(squeeze(met_images_coil_lim(:,:,slices,2,time_pts)), [0, max(met_images_coil_lim(:,:,slices,2,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (coil lim)');
imagescn(squeeze(met_images_coil_lim(:,:,slices,3,time_pts)), [0, max(met_images_coil_lim(:,:,slices,3,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

% multiresolution
slices = round(size(met_images_mres{1}, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (multires)');
imagescn(met_images_mres{1}(:,:,slices,time_pts), [0, max(met_images_mres{1}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (multires)');
imagescn(met_images_mres{2}(:,:,slices,time_pts), [0, max(met_images_mres{2}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (multires)');
imagescn(met_images_mres{3}(:,:,slices,time_pts), [0, max(met_images_mres{3}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

% noise
slices = round(size(met_images_mres_noise{1}, 3) / 2);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (noise)');
imagescn(met_images_mres_noise{1}(:,:,slices,time_pts), [0, max(met_images_mres_noise{1}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (noise)');
imagescn(met_images_mres_noise{2}(:,:,slices,time_pts), [0, max(met_images_mres_noise{2}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (noise)');
imagescn(met_images_mres_noise{3}(:,:,slices,time_pts), [0, max(met_images_mres_noise{3}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;



%% final met images (after mri)
slices = 1:size(met_images_upsampled{1}, 3);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (unified)');
imagescn(met_images_upsampled{1}(:,:,slices,time_pts), [0, max(met_images_upsampled{1}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (unified)');
imagescn(met_images_upsampled{2}(:,:,slices,time_pts), [0, max(met_images_upsampled{2}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (unified)');
imagescn(met_images_upsampled{3}(:,:,slices,time_pts), [0, max(met_images_upsampled{3}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

% AUCs
pyrAUC = sum(met_images_upsampled{1}, 4);
figure(Name='Pyr AUC (unified)');
imagescn(pyrAUC, [0 max(pyrAUC, [], 'all')], [1 numel(slices)]);
colormap fire;

lacAUC = sum(met_images_upsampled{2}, 4);
figure(Name='Lac AUC (unified)');
imagescn(lacAUC, [0 max(lacAUC, [], 'all')], [1 numel(slices)]);
colormap fire;

bicAUC = sum(met_images_upsampled{3}, 4);
figure(Name='Bic AUC (unified)');
imagescn(bicAUC, [0 max(bicAUC, [], 'all')], [1 numel(slices)]);
colormap fire;
