clear; close all;

% TISSUE STRUCTURE
mask = load('util/mask.mat').masks;
heart = tissue_structure("heart", mask, ["lv" "rv" "lvmy" "rvmy"]);
k_trans = [1, 1, 0.2, 0.4]; 
heart = heart.create_k_trans_map(k_trans);

% PK MODEL PARAMETERS
tr = 3.6;
n_t = 30;
flips = repmat([20; 30; 30;], 1, n_t) .* (pi/180);

n_compartments = 4;

% METABOLITE AND/OR TISSUE SPECIFIC PARAMETERS
mz0 = [0,1,0,0; % pyr
       0,0,0,0; % lac
       0,0,0,0]; % bic

r1 = [1/30, 1/25, 1/25]; % pyr, lac, bic

k = [0.0075, 0.0045, 0.06, 0.02; % lac
     0.0011, 0.0005, 0.04, 0.01]; % bic

t_arrival = [7, 0, 10, 14]; % lv, rv, lvmy, rvmy
t_bolus = 10;


% RUN PK MODEL
images = pk_model.run_pk_model(mz0, r1, k, flips, tr, heart, t_arrival=t_arrival, t_bolus=t_bolus);

% RUN MRI SYSTEM
sample_size = [32 32 11; 16 16 11; 24 24 11];
SNR = [150 40 20];
coil_lim = [0.4 1.2];
output_size = [32 32 11; 32 32 11; 32 32 11];
augmentation_params = struct(...
    "XTranslation", [-1,1], ...
    "YTranslation", [-1,1], ...
    "Scale", [0.95,1.1], ...
    "XReflection", true, ...
    "Rotation", [-5,5]);

met_images_mres = mri_system.run_mri_system(images, sample_size, SNR, coil_lim, heart.Mask, output_size, augmentation_params);

%% DISPLAY
figure;
imagescn(met_images_mres{1}(:,:,5,:), [0, max(met_images_mres{1}(:,:,5,:), [], 'all')])

figure;
imagescn(met_images_mres{2}(:,:,5,:), [0, max(met_images_mres{2}(:,:,5,:), [], 'all')])

figure;
imagescn(met_images_mres{3}(:,:,5,:), [0, max(met_images_mres{3}(:,:,5,:), [], 'all')])
