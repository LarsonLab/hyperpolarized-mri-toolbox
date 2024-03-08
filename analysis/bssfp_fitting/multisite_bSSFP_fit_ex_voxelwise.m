%% voxelwise fitting example

%% Rat Example: bSSFP-lactate Data 

clear; close all;

% load data
data_dir = './../../sample_data/Rat Kidneys Spiral bSSFP/';
load(strcat(data_dir,'shot1_bssfp_lac.mat'))
load(strcat(data_dir,'roi.mat'))

%pick roi to compute maps over
roi_name = kidneys;
% dilate mask before masking data to avoid kPL edge effects
roi_new = dilatemask3D(roi_name, ones([5 5]));

% mask bssfp data with roi, flatten and save indices
[nx, ny, nz, nt] = size(pyr_exp1);
roi_new_4d = repmat(roi_new, [1 1 1 nt]);
pyr_flat = reshape(pyr_exp1 .* roi_new_4d, [nx*ny*nz, nt]);
indices_bSSFP = find(pyr_flat(:,1));
lac_flat = reshape(lac_exp1 .* roi_new_4d, [nx*ny*nz, nt]);
pyr_flat = pyr_flat(any(pyr_flat ~= 0, 2), :); %remove pixels that are masked out
lac_flat = lac_flat(any(lac_flat ~= 0, 2), :); 
S = abs(cat(3, pyr_flat, lac_flat)); %use only pyr and lac for now
S = permute(S, [1 3 2]);

% ----set fitting params-----
% flip angles
FAP = 3; % [deg]
FAL = 60; % [deg]
Nz = 16; interleaves = 4;
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60]; % bSSFP catalyzation flip angles

% TR
TRP = .18; % [s]
TRL = .01529; % [s] 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz*interleaves]);
cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL]; % bSSFP catalyzation TR

% relaxation
params_fixed.R1P = 1/30; 
params_fixed.R1L = 1/25;
params_est.R2L = 1/0.8; 
params_est.kPL = 0.002;

TempRes = 4;
params_est.Mz0_P = 3;
params_est.Mz0_L = 0.5;
acq_sequence = ["3DGRE", "3DbSSFP"];
verbose = 0;

[fitparams, error, ~, ~] = multisite_bSSFP_fit(S, params_fixed, params_est, flips, TR, TempRes, acq_sequence,...
    'cat_flips', cat_flips, 'cat_TR', cat_TR, 'verbose', verbose);

%reshape kPLs and error maps back to 3D
kpl_map_flat = zeros([nx*ny*nz, 1]);
kpl_map_flat(indices_bSSFP) = fitparams.kPL';
kpl_map = reshape(kpl_map_flat, [nx ny nz]);

NRMSElac_map_flat = zeros([nx*ny*nz, 1]); % normalized root mean sq error for lac fit
NRMSElac_map_flat(indices_bSSFP) = error.NRMSE.L';
NRMSE_lac_map = reshape(NRMSElac_map_flat, [nx ny nz]);

if isfield(fitparams, 'R2L')
    T2L_map_flat = zeros([nx*ny*nz, 1]);
    T2L_map_flat(indices_bSSFP) = (1 ./ fitparams.R2L)';
    T2L_map = reshape(T2L_map_flat, [nx ny nz]);
end

% save maps if you'd like
% save_name = './../../../kpl_map.mat';
% if isfield(fitparams, 'R2L')
%     save(save_name, 'kpl_map', 'T2L_map', 'NRMSE_lac_map')
% else
%     save(save_name, 'kpl_map', 'NRMSE_lac_map')
% end


%% visualize maps on anatomical/localizer images
% modified from Shuyu's overlay_c13_proton.m

%load(save_name) %load saved maps

sl = 8:9; % slices to view
crop_idx = {65:165, 30:130}; % change this based on localizers
kpl_scale = [0 .01];
t2l_scale = [0 3];
nrmse_scale = [0 0.3];

%load localizers and roi
load(strcat(data_dir,'localizer.mat'))
load(strcat(data_dir,'roi.mat'))
mask = kidneys;

%plot kpl map
base_images = localizers(:,:,sl);
overlay_images = fermi_filter(kpl_map(:,:,sl),2); % can do the same with error maps
overlay_images = imresize(overlay_images,[size(base_images,1),size(base_images,2)],'nearest');
mask = imresize(mask,[size(base_images,1),size(base_images,2)],'nearest');
overlay_images = overlay_images .* mask(:,:,sl);
base_images2 = crop_images(base_images,crop_idx);
overlay_images = crop_images(overlay_images,crop_idx);
figure;imagescn_overlay(base_images2,[], overlay_images, kpl_scale, [1 numel(sl)], 1e-6, 0.4, 'fire', [], 1)

%plot T2L map
if isfield(fitparams, 'R2L')
    overlay_images = T2L_map(:,:,sl);
    overlay_images = imresize(overlay_images,[size(base_images,1),size(base_images,2)],'nearest');
    overlay_images = overlay_images .* mask(:,:,sl);
    overlay_images = crop_images(overlay_images,crop_idx);
    figure;imagescn_overlay(base_images2,[], overlay_images, t2l_scale, [1 numel(sl)], [], 0.3, 'fire', [], 1)
end

%plot NRMSE Lac map
overlay_images = NRMSE_lac_map(:,:,sl);
overlay_images = imresize(overlay_images,[size(base_images,1),size(base_images,2)],'nearest');
overlay_images = overlay_images .* mask(:,:,sl);
overlay_images = crop_images(overlay_images,crop_idx);
figure;imagescn_overlay(base_images2,[], overlay_images, nrmse_scale, [1 numel(sl)], [], 0.3, 'fire', [], 1)


%% util functions

function dilatedmask = dilatemask3D(mask, kernel)
% example use case: roi_new = dilatemask3D(roi_new, ones([5 5]));
% change size of ones matrix to control how dilated mask should be

    numsl = size(mask,3);
    dilatedmask = zeros(size(mask));
    for sl=1:numsl
        dilatedmask(:,:,sl) = imdilate(mask(:,:,sl),strel(kernel));
    end
end