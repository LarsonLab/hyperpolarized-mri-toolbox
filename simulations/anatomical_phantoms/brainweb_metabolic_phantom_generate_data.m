%% script used to generate metabolic phantom training data

%%
clear; close all;

%kinetic rates
kineticRates = [0,   randn_rng([.007 .035]); 
               0, randn_rng([0.002 0.007])]; % in order [vasc, GM, WM]
kineticRates(:,3) = kineticRates(:,2) + [randn_rng([-0.005 0.005]); 0];

%ktrans
ktransScales = [              1, randn_rng([.15 .25]);
                randn_rng([5 11]), randn_rng([.35 .45])];
diffWM = randn_rng([-0.005 -0.015]);
ktransScales(:,3) = ktransScales(:,2) + [diffWM; diffWM];

isFuzzy = true;
matSize = [32 32 8;
           16 16 8;
           16 16 8];
outputSize = [64 64 8];
linear_kTRANS_grad = true;
brain_idx = randi(19);

% define simulation parameters: Tarrival, Tbolus, TR, Nt, R1, flips,
% std_noise
simParams.Tarrival = randn_rng([-4 4], true);
simParams.Tbolus = 8;
simParams.TR = 4;
simParams.Nt = 20;
simParams.R1 = [1/30 1/25 1/25];
simParams.flips = repmat([20; 30; 30],[1 simParams.Nt])*pi/180;
simParams.SNR = [randn_rng([70 320]) randn_rng([15 75]) randn_rng([10 35])];
simParams.coil_lim = [randn_rng([0.2 0.6]) 1.2];

% define augmentation parameters
augmentParams.XTranslation = [-2 2];
augmentParams.YTranslation = [-2 2];
augmentParams.Scale = [0.95 1.2];
augmentParams.XReflection = true;
augmentParams.Rotation = [-5 5];

% input funciton and Mz0
input_function = realistic_input_function(simParams.Nt, simParams.TR, simParams.Tarrival, simParams.Tbolus);
Mz0P_scale = randn_rng([.5 1]); Mz0L_scale = randn_rng([0 .01]); Mz0B_scale = randn_rng([0 .005]);
Mz0 = [input_function(1), input_function(1)*Mz0P_scale,     input_function(1)*Mz0P_scale;
                       0, input_function(1)*Mz0L_scale,   input_function(1)*Mz0L_scale;
                       0, input_function(1)*Mz0B_scale,   input_function(1)*Mz0B_scale];

[k_trans, k_maps, Mz0_maps, metImages, w] = brainweb_metabolic_phantom(kineticRates, ktransScales, Mz0, matSize, outputSize, simParams, input_function, isFuzzy, linear_kTRANS_grad, augmentParams, brain_idx);

%% save to h5

savepath = "";
base_filename = "1_1";

for sl=1:size(k_trans,3)

    filepath = fullfile(savepath, strcat(base_filename,"_",num2str(sl),".h5"));

    k_trans_slice = squeeze(k_trans(:,:,sl));
    k_maps_slice = squeeze(k_maps(:,:,sl,:));
    Mz0_maps_slice = squeeze(Mz0_maps(:,:,sl,:));
    metImages_slice = squeeze(metImages(:,:,sl,:,:));
    weights_slice = squeeze(w(:,:,sl));
    
    metImages_slice = squeeze(permute(metImages_slice, [3 4 1 2]));
    
    % write data to h5
    h5create(filepath, "/kTRANS", size(k_trans_slice))
    h5write(filepath, "/kTRANS", k_trans_slice)
    h5create(filepath, "/kMaps", size(k_maps_slice))
    h5write(filepath, "/kMaps", k_maps_slice)
    h5create(filepath, "/Mz0Maps", size(Mz0_maps_slice))
    h5write(filepath, "/Mz0Maps", Mz0_maps_slice)
    h5create(filepath, "/metImages", size(metImages_slice))
    h5write(filepath, "/metImages", metImages_slice)
    h5create(filepath, "/weights", size(weights_slice))
    h5write(filepath, "/weights", weights_slice)
    
    % write attrs to h5
    h5writeatt(filepath,"/","kineticRates", kineticRates);
    h5writeatt(filepath,"/","ktransScales", ktransScales);
    h5writeatt(filepath,"/","SNR", simParams.SNR);
    h5writeatt(filepath,"/","Tarrival", simParams.Tarrival);
    h5writeatt(filepath,"/","coil_lim", simParams.coil_lim);
    h5writeatt(filepath,"/","Mz0_scale", [Mz0P_scale Mz0L_scale Mz0B_scale]);
    h5writeatt(filepath,"/","brain_idx", brain_idx);
end

%% VISUALIZATION all slices
% visualize kTRANS and k maps

slices = 1:size(k_trans,3);
tpts = 1:3:simParams.Nt;

% visualize kTRANS
figure,
imagescn(k_trans(:,:,slices),[0 max(k_trans(:,:,slices),[],'all')], [1 numel(slices)]); colormap fire;

% visualize kinetic rates
figure,
imagescn(k_maps(:,:,slices,1),[0 max(k_maps(:,:,slices,1),[],'all')], [1 numel(slices)]); colormap fire;

figure, 
imagescn(k_maps(:,:,slices,2),[0 max(k_maps(:,:,slices,2),[],'all')], [1 numel(slices)]); colormap fire;

% visualize Mz0 maps
figure,
imagescn(Mz0_maps(:,:,slices,1),[0 max(Mz0_maps(:,:,slices,1),[],'all')], [1 numel(slices)]); colormap fire;

figure, 
imagescn(Mz0_maps(:,:,slices,2),[0 max(Mz0_maps(:,:,slices,2),[],'all')], [1 numel(slices)]); colormap fire;

figure, 
imagescn(Mz0_maps(:,:,slices,3),[0 max(Mz0_maps(:,:,slices,3),[],'all')], [1 numel(slices)]); colormap fire;

% visualize metImages

%pyruvate
figure,
imagescn(squeeze(metImages(:,:,slices,1,tpts)),[0 max(squeeze(metImages(:,:,slices,1,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

%lactate
figure,
imagescn(squeeze(metImages(:,:,slices,2,tpts)),[0 max(squeeze(metImages(:,:,slices,2,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

%bicarb
figure,
imagescn(squeeze(metImages(:,:,slices,3,tpts)),[0 max(squeeze(metImages(:,:,slices,3,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

% visualize AUCs
%pyruvate
pyrAUC = sum(squeeze(metImages(:,:,slices,1,:)),length(size(squeeze(metImages(:,:,slices,1,:)))));
figure,
imagescn(pyrAUC,[0 max(pyrAUC,[],'all')], [1 numel(slices)]); colormap fire;

%lactate
lacAUC = sum(squeeze(metImages(:,:,slices,2,:)),length(size(squeeze(metImages(:,:,slices,2,:)))));
figure,
imagescn(lacAUC,[0 max(lacAUC,[],'all')], [1 numel(slices)]); colormap fire;

%bicarb
bicAUC = sum(squeeze(metImages(:,:,slices,3,:)),length(size(squeeze(metImages(:,:,slices,3,:)))));
figure,
imagescn(bicAUC,[0 max(bicAUC,[],'all')], [1 numel(slices)]); colormap fire;

% visualize voxel dynamics
sl = round(size(k_trans, 3)/2); %pick middle slice
metImages2 = squeeze(metImages(:,:,sl,:,:));

% visualize dynamics
S_1 = [squeeze(metImages2(18,30,1,:)) squeeze(metImages2(18,30,2,:)) squeeze(metImages2(18,30,3,:))];
S_2 = [squeeze(metImages2(38,24,1,:)) squeeze(metImages2(38,24,2,:)) squeeze(metImages2(38,24,3,:))];
S_3 = [squeeze(metImages2(37,50,1,:)) squeeze(metImages2(37,50,2,:)) squeeze(metImages2(37,50,3,:))];

figure('position', [0, 0, 300, 900]);
subplot(311)
plot(S_1, 'LineWidth', 2)
legend("pyr","lac","bic")
title("x=18, y=30")
subplot(312)
plot(S_2, 'Linewidth', 2)
title("x=38, y=24")
subplot(313)
plot(S_3, 'Linewidth', 2)
title("x=37, y=50")

%% VISUALIZATION one slice
% visualize kTRANS and k maps

slices = 4;
tpts = 1:3:simParams.Nt;

% visualize kTRANS and kinetic rates
%kTRANS
figure('position', [0, 0, 300, 900]);
subplot(311)
imagesc(squeeze(k_trans(:,:,slices)));
title("k_{TRANS}")
colormap fire;
axis off square;
colorbar;

%kPL
subplot(312)
imagesc(squeeze(k_maps(:,:,slices,1)));
title("k_{PL}")
colormap fire;
axis off square;
colorbar;

%kPB
subplot(313)
imagesc(squeeze(k_maps(:,:,slices,2)));
title("k_{PB}")
colormap fire;
axis off square;
colorbar;

% visualize Mz0 maps
%Mz0P
figure('position', [0, 0, 300, 900]);
subplot(311)
imagesc(squeeze(Mz0_maps(:,:,slices,1)));
title("Mz0P")
colormap fire;
axis off square;
colorbar;

%Mz0L
subplot(312)
imagesc(squeeze(Mz0_maps(:,:,slices,2)));
title("Mz0L")
colormap fire;
axis off square;
colorbar;

%Mz0B
subplot(313)
imagesc(squeeze(Mz0_maps(:,:,slices,3)));
title("Mz0B")
colormap fire;
axis off square;
colorbar;

% visualize metImages
%pyruvate
figure,
imagescn(squeeze(metImages(:,:,slices,1,tpts)),[0 max(squeeze(metImages(:,:,slices,1,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

%lactate
figure,
imagescn(squeeze(metImages(:,:,slices,2,tpts)),[0 max(squeeze(metImages(:,:,slices,2,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

%bicarb
figure,
imagescn(squeeze(metImages(:,:,slices,3,tpts)),[0 max(squeeze(metImages(:,:,slices,3,tpts)),[],'all')], [numel(slices) numel(tpts)]); colormap fire;

% visualize AUCs
%pyruvate
figure('position', [0, 0, 300, 900]);
subplot(311)
pyrAUC = sum(squeeze(metImages(:,:,slices,1,:)),length(size(squeeze(metImages(:,:,slices,1,:)))));
imagesc(pyrAUC);
title("Pyruvate AUC")
colormap fire;
axis off square;
colorbar;

subplot(312)
lacAUC = sum(squeeze(metImages(:,:,slices,2,:)),length(size(squeeze(metImages(:,:,slices,2,:)))));
imagesc(lacAUC);
title("Lactate AUC")
colormap fire;
axis off square;
colorbar;

subplot(313)
bicAUC = sum(squeeze(metImages(:,:,slices,3,:)),length(size(squeeze(metImages(:,:,slices,3,:)))));
imagesc(bicAUC);
title("Bicarbonate AUC")
colormap fire;
axis off square;
colorbar;

% visualize voxel dynamics
sl = round(size(k_trans, 3)/2); %pick middle slice
metImages2 = squeeze(metImages(:,:,sl,:,:));

% visualize dynamics
S_1 = [squeeze(metImages2(18,30,1,:)) squeeze(metImages2(18,30,2,:)) squeeze(metImages2(18,30,3,:))];
S_2 = [squeeze(metImages2(38,24,1,:)) squeeze(metImages2(38,24,2,:)) squeeze(metImages2(38,24,3,:))];
S_3 = [squeeze(metImages2(37,50,1,:)) squeeze(metImages2(37,50,2,:)) squeeze(metImages2(37,50,3,:))];

figure('position', [0, 0, 300, 900]);
subplot(311)
plot(S_1, 'LineWidth', 2)
legend("pyr","lac","bic")
title("x=18, y=30")
subplot(312)
plot(S_2, 'Linewidth', 2)
title("x=38, y=24")
subplot(313)
plot(S_3, 'Linewidth', 2)
title("x=37, y=50")


