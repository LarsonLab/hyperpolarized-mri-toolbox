clear; close all;
addpath('../');
addpath('../../pk_models/');
addpath('../../../utilities/')

%% PARAMETERS ----------------------------------------------------
% tissue
mask = load('../util/cardiac/cardiac_mask.mat').masks;

k_trans = [1, 1, 0.2, 0.4]; 

% pk model params
tr = 3;
n_t = 20;
flips = repmat([20; 30; 30], 1, n_t) .* (pi/180); 
r1 = [1/30 1/25 1/25];
k = [0.0075, 0.0045, 0.06, 0.02;
     0.0011, 0.0005, 0.0400, 0.01]; % in order [lv, rv, lv_mc, rv_mc].

t_arrival = [7,0,10,14];
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

mz0 = repmat(reshape(input_function(:, 1), 1, n_tissues), 3, 1) .* mz0_constants; % TODO: check this. it's a little convoluted and lots of magic numbers. also, lowk a bad way to do this anyways

% mri
coil_lim = [0.4 1.2];
sample_size = [32,32,5; 16,16,5; 16,16,5];
SNR = [150 40 20];
output_size = [32,32,5];

%% RUNNING THE MODEL ----------------------------------------------
% tissue
heart = tissue_structure("heart", mask, ["lv" "rv" "lvmy" "rvmy"]);
heart = heart.create_k_trans_map(k_trans);

% pk model
images = pk_model.run_pk_model(mz0, r1, k, flips, tr, heart, input_function=input_function);

% mri
met_images_mres = mri_system.run_mri_system(images, sample_size, SNR, coil_lim, heart.Mask, output_size);

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

% AUCs
pyrAUC = sum(met_images_mres{1}, 4);
figure(Name='Pyr AUC (unified)');
imagescn(pyrAUC, [0 max(pyrAUC, [], 'all')], [1 numel(slices)]);
colormap hot;

lacAUC = sum(met_images_mres{2}, 4);
figure(Name='Lac AUC (unified)');
imagescn(lacAUC, [0 max(lacAUC, [], 'all')], [1 numel(slices)]);
colormap hot;

bicAUC = sum(met_images_mres{3}, 4);
figure(Name='Bic AUC (unified)');
imagescn(bicAUC, [0 max(bicAUC, [], 'all')], [1 numel(slices)]);
colormap hot;
