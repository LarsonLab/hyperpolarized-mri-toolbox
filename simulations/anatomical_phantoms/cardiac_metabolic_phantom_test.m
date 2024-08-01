% Quick testing script for the cardiac metabolic phantom 
clear; close all;

kineticRates = [0.0075, 0.0045, 0.0264, 0.0179;
                0.0011, 0.0005, 0.0100, 0.0017]; % in order [lv, rv, lv_mc, rv_mc].
ktransScales = [1.2, 1, 0.2, 0.2]; 

matSize = [32,32,5];

heart_idx = 1;

% define simulation parameters: Tarrival, Tbolus, TR, Nt, R1, flips,
% std_noise
simParams.Tarrival = [5,0,10,10]; % in order lv, rv, lvmy, rvmy
simParams.Tbolus = 1;
simParams.TR = 3; % changes over time, but this seems like a pretty good estimate
% TR = 3 * 60/heart_rates(I_subject);
% subject_ids = [6  7  8  9  10 11 13 18];
% heart_rates = [76 80 69 64 50 77 61 57];
simParams.Nt = 20;
simParams.R1 = [1/30 1/25 1/25];
simParams.flips = repmat([20; 30; 30],[1 simParams.Nt])*pi/180; 
simParams.SNR = [150 40 20];
simParams.coil_lim = [0.4 1.2]; % TODO: account for coil sensitivity

Mz0_constants = [1, 1, 0.5, 0.5;
                0, 0, 0.01, 0.01;
                0, 0, 0.005, 0.005];

[k_trans, k_maps, Mz0_maps, metImages] = cardiac_metabolic_phantom(kineticRates, ktransScales, Mz0_constants, matSize, simParams, heart_idx);

slices = 1:size(k_trans,3);
tpts = 1:simParams.Nt;

%% visualize kTRANS and k maps

% visualize kTRANS
figure("Name","kTRANS")
imagescn(k_trans(:,:,slices),[0 max(k_trans(:,:,slices),[],'all')], [1 numel(slices)]); colormap hot;

% visualize kinetic rates
figure("Name","kPL");
imagescn(k_maps(:,:,slices,1),[0 max(k_maps(:,:,slices,1),[],'all')], [1 numel(slices)]); colormap hot;

figure("Name","kPB");
imagescn(k_maps(:,:,slices,2),[0 max(k_maps(:,:,slices,2),[],'all')], [1 numel(slices)]); colormap hot;


%% visualize Mz0 maps
figure("Name","Mz0 maps Pyruvate");
scale = [0 max(Mz0_maps(:,:,slices,1),[],'all')];
if scale(2) == 0
    scale(2) = 1;
end
imagescn(Mz0_maps(:,:,slices,1), scale, [1 numel(slices)]); colormap hot;

figure("Name","Mz0 maps Lactate");
scale = [0 max(Mz0_maps(:,:,slices,2),[],'all')];
if scale(2) == 0
    scale(2) = 1;
end
imagescn(Mz0_maps(:,:,slices,3), scale, [1 numel(slices)]); colormap hot;


figure("Name","Mz0 maps Bicarb");
scale = [0 max(Mz0_maps(:,:,slices,3),[],'all')];
if scale(2) == 0
    scale(2) = 1;
end
imagescn(Mz0_maps(:,:,slices,3), scale, [1 numel(slices)]); colormap hot;


%% visualize metImages
%pyruvate
figure("Name","Pyruvate Met Images");
imagescn(squeeze(metImages(:,:,slices,1,tpts)),[0 max(squeeze(metImages(:,:,slices,1,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap hot;

%lactate
figure("Name","Lactate Met Images");
imagescn(squeeze(metImages(:,:,slices,2,tpts)),[0 max(squeeze(metImages(:,:,slices,2,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap hot;

%bicarb
figure("Name","Bicarb Met Images");
imagescn(squeeze(metImages(:,:,slices,3,tpts)),[0 max(squeeze(metImages(:,:,slices,3,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap hot;


%% visualize AUCs

%pyruvate
pyrAUC = sum(squeeze(metImages(:,:,:,1,:)),length(size(squeeze(metImages(:,:,:,1,:)))));
figure("Name","Pyruvate AUC");
imagescn(pyrAUC,[0 max(pyrAUC,[],'all')], [1 numel(slices)]); colormap hot;


%lactate
lacAUC = sum(squeeze(metImages(:,:,:,2,:)),length(size(squeeze(metImages(:,:,:,2,:)))));
figure("Name","Lactate AUC");
imagescn(lacAUC,[0 max(lacAUC,[],'all')], [1 numel(slices)]); colormap hot;


%bicarb
bicAUC = sum(squeeze(metImages(:,:,:,3,:)),length(size(squeeze(metImages(:,:,:,3,:)))));
figure("Name","Bicarb AUC");
imagescn(bicAUC,[0 max(bicAUC,[],'all')], [1 numel(slices)]); colormap hot;
