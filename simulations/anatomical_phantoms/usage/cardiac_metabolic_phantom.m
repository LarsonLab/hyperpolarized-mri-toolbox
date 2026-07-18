clear; close all;
addpath('../'); % phantom functions
addpath('../../pk_models/'); % realistic_input_function, simulate_Nsite_model
addpath('../../../utilities/'); % fire

%% PARAMETERS ----------------------------------------------------
% tissue
mask = load('../util/cardiac/cardiac_mask_10.mat').masks;

k_trans = [1, 1, 0.2, 0.4]; 

% pk model params
tr = 3.6;
n_t = 30;
flips = repmat([20; 30; 30], 1, n_t) .* (pi/180); 
r1 = [1/30 1/25 1/20];
k = [0.013, 0.010, 0.025, 0.020;
     0.001, 0.001, 0.010, 0.001]; % in order [lv, rv, lv_mc, rv_mc].

t_arrival = [6,0,8,2];
t_bolus = 1;

% input function and mz0
n_tissues = size(k, 2);
n_mets = size(flips, 1);

input_function = zeros(n_tissues, n_t);
for i_tissue = 1:n_tissues
    input_function(i_tissue, :) = realistic_input_function(n_t, tr, t_arrival(i_tissue), t_bolus);
end

mz0_constants = [1, 1, 0.5, 0.5;
                0, 0, 0.01, 0.01;
                0, 0, 0.005, 0.005];

mz0 = repmat(reshape(input_function(:, 1), 1, n_tissues), 3, 1) .* mz0_constants;

% mri
coil_lim = [0.4 1.2];
sample_size = [25,25,5; 13,13,5; 13,13,5];
snr = [220 70 12];
output_size = [32,32,5];

%% RUNNING THE MODEL ----------------------------------------------
% tissue
heart = tissue_structure("heart", mask, ["lv" "rv" "lvmy" "rvmy"]);
heart = heart.create_k_trans_map(k_trans);

% pk model
images = pk_model.run_pk_model(mz0, r1, k, flips, tr, heart, input_function=input_function);

% mri
met_images_mres = mri_system.run_mri_system(images, sample_size, snr, coil_lim, heart, output_size, include_bg_noise=false);

%% DISPLAY --------------------------------------------------------
slices = 1:5:46;
figure(Name='kTRANS');
imagescn(heart.K_trans_map(:,:,slices), [0 max(heart.K_trans_map(:,:,slices), [], 'all')], [1 numel(slices)]);
colormap hot;

slices = 1:size(met_images_mres{1}, 3);
time_pts = 1:n_t;
figure(Name='Pyruvate');
imagescn(met_images_mres{1}(:,:,slices,time_pts), [0, max(met_images_mres{1}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap hot;

figure(Name='Lactate');
imagescn(met_images_mres{2}(:,:,slices,time_pts), [0, max(met_images_mres{2}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap hot;

figure(Name='Bicarb (unified)');
imagescn(met_images_mres{3}(:,:,slices,time_pts), [0, max(met_images_mres{3}(:,:,slices,time_pts), [], 'all')], [numel(slices) numel(time_pts)]);
colormap hot;

%% AUCs
pyrAUC = sum(met_images_mres{1}, 4);
lacAUC = sum(met_images_mres{2}, 4);
bicAUC = sum(met_images_mres{3}, 4);

% AUC ratios
lac_to_pyr = lacAUC ./ pyrAUC;
bic_to_pyr = bicAUC ./ pyrAUC;

% remove the background noise
parts_to_keep = sum(heart.Mask(), 4);
parts_to_keep = imresize3(parts_to_keep, output_size, "nearest");
lac_to_pyr(parts_to_keep == 0) = 0;
bic_to_pyr(parts_to_keep == 0) = 0;

figure(Name='Lac/Pyr AUC (unified)');
imagescn(lac_to_pyr, [0, 0.5], [1 numel(slices)]);
colormap hot;

figure(Name='Bic/Pyr AUC (unified)');
imagescn(bic_to_pyr, [0, 0.15], [1 numel(slices)]);
colormap hot;


