clear; close all;
addpath('../');
addpath('../../pk_models/');
addpath('../../../utilities/')

%% PARAMETERS
% tissue
mask = load('../util/brainweb/1/brainweb_fuzzy.mat').im_mask;

k_trans = [1, 0.2, 0.2;
           3, 0.4, 0.4];

% pk model params
tr = 3;
n_t = 20;
flips = repmat([20; 30; 30;], 1, n_t) .* (pi/180);
r1 = [1/30, 1/25, 1/20]; % pyr, lac, bic
k = [0, 0.03, 0.025; % lac, in order [vasc, gm, wm]
     0, 0.01, 0.005]; % bic
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
coil_lim = [0.2, 0.6];
augment_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5], ...
    "ZTranslation", [-20, 20]);

sample_size = [32 32 8; 16 16 8; 16 16 8];
snr = [150 40 20];
output_size = [64 64 8];

%% RUNNING THE MODEL -----------------------------------------------------------

% tissue
brain = tissue_structure("brain", mask, ["vasc", "gm", "wm"], [100,100,100]);
brain = brain.create_k_trans_map(k_trans);

% pk model
images = pk_model.run_pk_model(mz0, r1, k, flips, tr, brain, input_function=input_function);

% mri
met_images_mres = mri_system.run_mri_system(images, sample_size, snr, coil_lim, brain, output_size, augment_params);

%% DISPLAY ---------------------------------------------------------------------
slices = 10:5:40;
figure(Name='kTRANS');
imagescn(brain.K_trans_map(:,:,slices), [0 max(brain.K_trans_map(:,:,slices), [], 'all')], [1 numel(slices)]);
colormap fire;

k_trans_dwnszd = imresize3(imresize3(brain.K_trans_map, [16 16 8]), [32 32 8]);
slices = 1:size(k_trans_dwnszd, 3);
figure(Name='kTRANS downsized');
imagescn(k_trans_dwnszd(:,:,slices), [0 max(k_trans_dwnszd(:,:,slices), [], 'all')], [1 numel(slices)]);
colormap fire;
%%
slices = 1:size(met_images_mres{1}, 3);
time_pts = 1:3:n_t;
figure(Name='Pyruvate (unified)');
imagescn(met_images_mres{1}(:,:,slices,time_pts), [0, max(met_images_mres{1}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Lactate (unified)');
imagescn(met_images_mres{2}(:,:,slices,time_pts), [0, max(met_images_mres{2}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;

figure(Name='Bicarb (unified)');
imagescn(met_images_mres{3}(:,:,slices,time_pts), [0, max(met_images_mres{3}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap fire;
%%
% AUCs
pyrAUC = sum(met_images_mres{1}, 4);
figure(Name='Pyr AUC (unified)');
imagescn(pyrAUC, [0 max(pyrAUC, [], 'all')], [1 numel(slices)]);
colormap fire;

lacAUC = sum(met_images_mres{2}, 4);
figure(Name='Lac AUC (unified)');
imagescn(lacAUC, [0 max(lacAUC, [], 'all')], [1 numel(slices)]);
colormap fire;

bicAUC = sum(met_images_mres{3}, 4);
figure(Name='Bic AUC (unified)');
imagescn(bicAUC, [0 max(bicAUC, [], 'all')], [1 numel(slices)]);
colormap fire;
