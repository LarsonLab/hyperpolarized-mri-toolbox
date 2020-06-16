% *** Tensor Image Enhancement and Optimal Array Combination ***
% Demos the use of high-order SVD (tensors) to enhance SNR of hyperpolarized 
% 13C MRSI data and optimal combination of array receivers.
% +++ Note: Peak index/width needs to be manually specified in this main 
% code(spec_idx_met) and 2d_epsi_c13.peak +++
% Requires MATLAB R2015+ and support of unix command environment

% References:
% [1] Brender JR, Kishimoto S, Merkle H, et al. Dynamic Imaging of Glucose  
%     and Lactate Metabolism by 13 C-MRS without Hyperpolarization. Sci Rep 
%     2019; 9(1):1-14.
% [2] Chen H-Y, Autry AW, Brender JR, et al. Tensor Image Enhancement and 
%     Optimal Multichannel Receiver Combination Analyses For Human Hyper-
%     polarized 13C MRSI. Magn Reson Med 2020; doi: 10.1002/mrm.28328
% [3] Vareth M, Lupo J, Larson P, et al. A comparison of coil combination
%     strategies in 3D multichannel MRSI reconstruction for patients with 
%     brain tumors. NMR Biomed 2018;31(11):e3929
%
% Licensed under CC BY-NC-SA 4.0.

%%
clear all
close all
clc

% assuming input is reconstructed HP 13C-dynamic MRSI data
addpath(genpath('utility'));
% set parameters
spec_idx_met = {180:250,1:50,110:145}; % spectral range of {pyr,lac,bic} in points
autoRank_on = true; % turn on information-based optimal rank selection (time-consuming)
output_Directory = 'outputSpectrum_TRI';


%% process pyruvate lactate
% **link to the sample dataset** https://zenodo.org/record/3895204
load('data_2DEPSI_dyno.mat');
% ---- noise prewhitening - prep step for channel combination ----
% generate noise decorrelation matrix
dmtx = generate_dmtx_recondata(csi_phased);
% prewhitening
csi_phased = noise_prewhitening(csi_phased, dmtx); 
Data_5D_orig = squeeze(permute(csi_phased,[4 1 3 2 5 6]));
size_Data_5D_orig = size(Data_5D_orig); % [f,x,y,coil,dyn]

% ---- B0 correction ----
peaksearch_range = spec_idx_met{1}; % use location of pyruvate peak to align B0 shifts
[Data_5D_orig,circ_shift_dist]= pb_correct_chemshift_brain(Data_5D_orig,peaksearch_range);

% ---- SVD coil combination ----
% Estimate Sensitivities
svdWeights = pb_est_coil_sensitivity(Data_5D_orig);

% Ideal channel combination using SVD coefficients
Data_5D_orig = Data_5D_orig .* ...
    repmat(permute(svdWeights,[4 2 3 1 5]),[size_Data_5D_orig(1) 1 1 1 size_Data_5D_orig(5)]);
Data_5D_orig = squeeze(sum(Data_5D_orig,4));

% % ---- Tensor Image Enhancement ----
% Determine Optimal Rank
if autoRank_on
    % try resample data
    sampOrig = 59;
    sampRecon = 256;
    BVratio = sampRecon/sampOrig;
    spec_idx_met_resamp = {41:58,1:12};
    Data_5D_resamp = reshape(resample(Data_5D_orig,sampOrig,sampRecon),...
        [sampOrig size_Data_5D_orig(2:3) size_Data_5D_orig(5)]);
    [sizeLR, CF_Params, Td_Params] = pb_autorank(Data_5D_resamp,spec_idx_met_resamp,BVratio,1);
    save('autoRank_params.mat','sizeLR','CF_Params','Td_Params');
else
    sizeLR = [6, 8, 10, 10];
end
%% Tensor thresholding
Data_5D_trunc = td_regularization(Data_5D_orig,sizeLR,0);

% ---- phase and baseline correction ----
Data_5D_pl_pb = pb_phase_baseline(Data_5D_orig,spec_idx_met);

%% process bicarbonate
% note aliasing of bicarb in EPSI spectrum
filepath = 'data_2DEPSI_dyno_bicarb.mat';
bicarb_params.spec_idx_met = {188:234,9:41,90:117}; % note new bicarb freq due to demodulation
bicarb_params.dmtx = dmtx;
bicarb_params.circ_shift_dist = circ_shift_dist;
bicarb_params.svdWeights = svdWeights;
sizeLR_bicarb = sizeLR;
met_specAUC_bicarb = pb_process_bicarb(filepath,bicarb_params,Data_5D_pl_pb,sizeLR_bicarb);

%% quantify metabolites
met_specAUC_before = spec2img(Data_5D_orig,spec_idx_met);
met_specAUC_after = spec2img(Data_5D_trunc,spec_idx_met);
met_specAUC_before{3} = met_specAUC_bicarb{1};
met_specAUC_after{3} = met_specAUC_bicarb{2};
pb_plot_temporal(met_specAUC_before,met_specAUC_after);

%% save processed spectrum
% note saved spectrum is demodulated at pyruvate/lactate frequency
saveRawSpec = true;
met_specAUC_raw = met_specAUC_before;
met_specAUC_TRI = met_specAUC_after;
Data_5D_TRI = Data_5D_trunc;
mkdir(output_Directory)
if saveRawSpec
    save(sprintf('%s/outputSpectrum_TRI.mat',output_Directory),'Data_5D_TRI','met_specAUC_raw','met_specAUC_TRI');
else
    save(sprintf('%s/outputSpectrum_TRI.mat',output_Directory),'Data_5D_TRI','met_specAUC_TRI');
end
