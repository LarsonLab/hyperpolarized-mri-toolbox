clear; close all;
addpath('../');

%% PARAMETERS
% tissue
mask = load('brainweb.mat').im_mask;

newmask = zeros(100, 100, 50, 3);
for i = 1:3
    newmask(:,:,:,i) = imresize3(mask(:,:,:,i), [100,100,50]);
end
mask = newmask;

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
n_compartments = 3;
input_function = zeros(n_compartments, n_t);
for i_cmp = 1:n_compartments
    input_function(i_cmp, :) = realistic_input_function(n_t, tr, t_arrival(i_cmp), t_bolus);
end
mz0 = [input_function(1), input_function(1)*.5, input_function(1)*.5;
       0, input_function(1)*.01, input_function(1)*.01;
       0, input_function(1)*.005, input_function(1)*.005];

% mri
coil_lim = [0.4, 1.2];
augmentation_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5]);

sample_size = [16 16 16; 16 16 16; 8 8 8];
SNR = [150 40 20];
output_size = [32 32 8; 32 32 8; 32 32 8];

%% RUNNING THE MODEL -----------------------------------------------------------

% tissue
brain = tissue_structure("brain", mask, ["vasc", "gm", "wm"]);
brain = brain.create_k_trans_map(k_trans);

% pk model
images = pk_model.run_pk_model(mz0, r1, k, flips, tr, brain, input_function=input_function);

% mri
met_images_mres = mri_system.run_mri_system(images, sample_size, SNR, coil_lim, brain.Mask, output_size, augmentation_params);

%% DISPLAY ---------------------------------------------------------------------
figure;
imagescn(met_images_mres{1}(:,:,5,:), [0, max(met_images_mres{1}(:,:,5,:), [], 'all')]);
colormap fire;

figure;
imagescn(met_images_mres{2}(:,:,5,:), [0, max(met_images_mres{2}(:,:,5,:), [], 'all')]);
colormap fire;

figure;
imagescn(met_images_mres{3}(:,:,5,:), [0, max(met_images_mres{3}(:,:,5,:), [], 'all')]);
colormap fire;
