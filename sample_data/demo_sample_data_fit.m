
clear all

%% Load data


I_sample_data = 1;

switch I_sample_data
    case 1
        datafile = 'Rat Kidneys EPI/exp1_vfa/exp1_vfa.mat'; zplot = 1; SNR_thresh = 1e5;
        phase_flag = 0;
    case 2
        datafile = 'Rat Kidneys EPI/exp2_constant/exp2_constant.mat'; zplot = 1; SNR_thresh = 1e5;
        phase_flag = 0;
    case 3
        datafile = 'TRAMP Mouse multislice EPI/exp_tramp_epi.mat'; zplot = 8; SNR_thresh = 1e5;
        phase_flag = 1;
end

load(datafile)
flips_all = [flips_pyr(:).';flips_lac(:).'] *pi/180;
size_pyr= size(pyr);
Nx_dims = length(size_pyr)-1;
% convert to [x,y,z,t]
if Nx_dims ==2
    pyr = reshape(pyr, [size_pyr(1:Nx_dims), 1, size_pyr(end)]);
    lac = reshape(lac, [size_pyr(1:Nx_dims), 1, size_pyr(end)]);
end

size_pyr= size(pyr);

% phase data based on dynamics
if phase_flag
    for x = 1:size_pyr(1)
        for y = 1:size_pyr(2)
            for z= 1:size_pyr(3)
                pyr(x,y,z,:) =  pyr(x,y,z,:) * exp(i*find_phase_corr(squeeze(pyr(x,y,z,:))));
                lac(x,y,z,:) =  lac(x,y,z,:) * exp(i*find_phase_corr(squeeze(lac(x,y,z,:))));
            end
        end
    end
    
end

AUC_pyr = real(sum(pyr,4));
AUC_lac = real(sum(lac,4));
AUC_ratio = AUC_lac./AUC_pyr;

data = cat(4, reshape(real(pyr), [size_pyr(1:3), 1, size_pyr(4)]), ...
    reshape(real(lac), [size_pyr(1:3), 1, size_pyr(4)]) );




%% Fit data

clear params_fixed params_est params_fit
R1P_est = 1/30; R1L_est = 1/25;  kPL_est = 0.01;
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

SNRmask = AUC_pyr > SNR_thresh;
It_fit = find(flips_lac > 1);

[params_fit Sfit] = fit_kPL(double(data(:,:,:,:,It_fit)), TR, flips_all(:, It_fit), params_fixed, params_est);

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
