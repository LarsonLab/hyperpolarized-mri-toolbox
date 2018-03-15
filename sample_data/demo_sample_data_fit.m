
clear all

%% Load data

% uncomment data to test
%datafile = 'TRAMP Mouse multislice EPI/exp_tramp_epi.mat'; zplot = 10; SNR_thresh = 1e5;
datafile = 'Rat Kidneys EPI/exp1_vfa/exp1_vfa.mat'; zplot = 1; SNR_thresh = 1e5;
%datafile = 'Rat Kidneys EPI/exp2_constant/exp2_constant.mat'; zplot = 1; SNR_thresh = 1e5;

load(datafile)
flips_all = [flips_pyr(:).';flips_lac(:).'] *pi/180;
size_pyr= size(pyr);
Nx_dims = length(size_pyr)-1;

AUC_pyr = sum(pyr,Nx_dims+1);
AUC_lac = sum(lac,Nx_dims+1);
AUC_ratio = AUC_lac./AUC_pyr;

pyr = reshape(pyr, [size_pyr(1:Nx_dims), 1, size_pyr(end)]);
lac = reshape(lac, [size_pyr(1:Nx_dims), 1, size_pyr(end)]);
data = cat(Nx_dims+1, pyr, lac);




%% Fit data

clear params_fixed params_est params_fit
R1P_est = 1/30; R1L_est = 1/25;  kPL_est = 0.01;
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

SNRmask = AUC_pyr > SNR_thresh;
It_fit = find(flips_lac > 1);

[params_fit Sfit] = fit_kPL(double(data(:,:,:,It_fit)), TR, flips_all(:, It_fit), params_fixed, params_est);

figure
subplot(221)
imagesc(AUC_pyr(:,:,zplot))
title('AUC pyruvate')
subplot(222)
imagesc(AUC_lac(:,:,zplot))
title('AUC lactate')
subplot(223)
imagesc(AUC_ratio(:,:,zplot) .* SNRmask(:,:,zplot)), colorbar
title('AUC lactate / AUC pyruvate')
subplot(224)
imagesc(params_fit.kPL(:,:,zplot) .* SNRmask(:,:,zplot)), colorbar
title('Fit k_{PL}')
