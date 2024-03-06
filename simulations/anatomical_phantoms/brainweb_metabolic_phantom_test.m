% Quick testing script for the brain web basic metabolic phantom 
clear; close all;

kineticRates = [0, 0.05, 0.03; 
                0,    0,    0]; % no kPB, effectively, in order [vasc, GM, WM]
ktransScales = [1, 0.3 0.3];
isFuzzy = false;
matSize = [16 16 8];

% define simulation parameters: Mz0, Tarrival, Tbolus, TR, Nt, R1, flips,
% std_noise
simParams.Mz0 = [0 0 0];
simParams.Tarrival = 0;
simParams.Tbolus = 8;
simParams.TR = 4;
simParams.Nt = 30;
simParams.R1 = [1/30 1/25 1/25];
simParams.flips = ones(3,simParams.Nt)*20*pi/180;
simParams.std_noise = 0.003; 

[k_trans, k_maps, metImages] = brainweb_metabolic_phantom(kineticRates, ktransScales, isFuzzy, matSize, simParams);

%% visualize kTRANS and k maps

slices = 1:matSize(3);
tpts = 1:3:simParams.Nt;

% visualize kTRANS
figure,
imagescn(k_trans(:,:,slices),[0 max(k_trans(:,:,slices),[],'all')], [1 numel(slices)]); colormap fire;

% visualize kinetic rates
figure,
imagescn(k_maps(:,:,slices,1),[0 max(k_maps(:,:,slices,1),[],'all')], [1 numel(slices)]); colormap fire;

%figure, % picking an arbitrarily small number for the lower scale bound since k13 == 0
%imagescn(k_maps(:,:,slices,2),[0 max(k_maps(:,:,slices,2),[],'all')], [1 numel(slices)]); colormap fire;

%% visualize metImages

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
pyrAUC = sum(squeeze(metImages(:,:,:,1,:)),4);
figure,
imagescn(pyrAUC,[0 max(pyrAUC,[],'all')], [1 numel(slices)]); colormap fire;

%lactate
lacAUC = sum(squeeze(metImages(:,:,:,2,:)),4);
figure,
imagescn(lacAUC,[0 max(lacAUC,[],'all')], [1 numel(slices)]); colormap fire;

%bicarb
bicAUC = sum(squeeze(metImages(:,:,:,3,:)),4);
figure,
imagescn(bicAUC,[0 max(bicAUC,[],'all')], [1 numel(slices)]); colormap fire;
