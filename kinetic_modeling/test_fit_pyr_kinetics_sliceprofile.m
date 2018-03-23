% Script for testing fit_kPL kinetic model fitting functions

clear all

% Test values
Tacq = 48; Nflips = 8;
TR = 3; TRall = TR/Nflips; Nall = Tacq/TRall; Nt = Tacq/TR;
R1P = 1/25; R1L = 1/25; R1B = 1/15; R1A = 1/25;
kPL = 0.05; kPB = 0.03; kPA = 0.02;
std_noise = 0.01;

% gamma variate input function - most realistic
tall = [1:Nall]*TRall;
Tarrival = 0; Tbolus = 12;
A = 4; B = Tbolus/4;
input_function = gampdf(tall-Tarrival,A,B);  % gamma distribution -continuous input
input_function = input_function/sum(input_function);% normalize for a total magnetization input = 1
Mz0 = zeros(4,1);

load('correction_factors_profiles.mat', 'zloc', 'flip_profile')
lo = find(zloc > -20, 1);
hi = find(zloc > 20, 1);
slice_profile = flip_profile(lo:5:hi); 
%slice_profile = 1; % no slice profile correction
M = length(slice_profile);

% Test over multiple combinations of flip angle schemes
k12 = 0.03; % for variable flip angle designs
flips(1:2,1:Nall,1) = ones(2,Nall)*30*pi/180  / sqrt(Nflips);  % constant, single-band
flips(1:2,1:Nall,2) = repmat([20;35]*pi/180,[1 Nall]) / sqrt(Nflips);  % constant, multi-band
flips(1:2,1:Nall,3) = [vfa_const_amp(Nall, pi/2, exp(-TRall * ( k12))); ... % max lactate SNR variable flip angle
    vfa_opt_signal(Nall, exp(-TRall * ( R1L)))];

flips = cat(1, flips, repmat( flips(2,:,:), [2 1 1]));

N_flip_schemes = size(flips,3);

figure
subplot(121) , plot(tall, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(tall, squeeze(flips(2,:,:))*180/pi)
title('Product (Lactate, Alanine, Bicarbonate) flips')
legend('constant','multiband',  'max lactate SNR vfa'); %, 'vfa', 'T1-effective vfa', 'Saturation Recovery')

% generate simulated data
noise_S = randn([4 Nall])*std_noise;  % same noise for all flip schedules
for Iflips = 1:N_flip_schemes
	for m = 1:M
        [Mxytemp(1:4, 1:Nall, m), Mz] = simulate_Nsite_model(Mz0, [R1P R1L R1B R1A], [kPL 0; kPB 0; kPA 0], flips(:,:,Iflips)*slice_profile(m) , TRall, input_function);
    end
    
    if m > 1
        % average over slice profile
        Mxy(1:4, 1:Nall, Iflips) = mean(Mxytemp,3);
    else
        Mxy(1:4, 1:Nall, Iflips) = Mxytemp;
    end
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy(:,:,Iflips) + noise_S;
end

    % average across Nflips per image ( multiple phase encodes )
    Mxy = squeeze( sum(reshape(Mxy, [4, Nflips, Nt, N_flip_schemes]),2) ) ;
    Sn = squeeze( sum(reshape(Sn, [4, Nflips, Nt, N_flip_schemes]),2) ) ;
std_noise = std_noise * sqrt(Nflips);

% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; R1B_est = 1/15; R1A_est = 1/25;
kPL_est = .02; kPB_est = .02; kPA_est = .02;
plot_fits = 1;

%% Test fitting
disp('Fitting kPL, kPB, and kPA, with fixed relaxation rates:')
disp('Including slice profile corrections')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], slice_profile, plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], slice_profile, plot_fits);
    
    % magnitude fitting with noise
   [params_fitn_mag(:,Iflips) Snfit_mag(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, slice_profile, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])
disp(['KPB  = ' num2str(kPB)])
disp(['KPA  = ' num2str(kPA)])

Iskip=9;
disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
k_fit = k_fit([[1:3],[1:3]+Iskip,[1:3]+2*Iskip]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:3],[1:3]+Iskip,[1:3]+2*Iskip]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:3],[1:3]+Iskip,[1:3]+2*Iskip]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])

t = [0:Nt-1]*TR ;

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits')% (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')
subplot(224) , plot(t, squeeze(Sn(4,:,:)))
hold on, plot(t, squeeze(Snfit_complex(3,:,:)),':')
plot(t, squeeze(Snfit_mag(3,:,:)),'--'), hold off
title('Alanine signals and fits')

return
%% Test fitting without slice profile correction
disp('Fitting kPL, kPB, and kPA, with fixed relaxation rates:')
disp('WITHOUT slice profile corrections')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], [], plot_fits);
    
    % magnitude fitting with noise
 %   [params_fitn_mag(:,Iflips) Snfit_mag(1:3,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])
disp(['KPB  = ' num2str(kPB)])
disp(['KPA  = ' num2str(kPA)])

Iskip=9;
disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
k_fit = k_fit([[1:3],[1:3]+Iskip,[1:3]+2*Iskip]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:3],[1:3]+Iskip,[1:3]+2*Iskip]); 
disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])
% disp('Noisy magnitude fit results:')
% k_fit = struct2array(params_fitn_mag);
% k_fit = k_fit([[1:3],[9:11],[17:19]]); 
% disp([['KPL  = ';'KPB  = ';'KPA  = '],num2str(reshape(k_fit,[3 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
%plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits')% (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
%plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')
subplot(224) , plot(t, squeeze(Sn(4,:,:)))
hold on, plot(t, squeeze(Snfit_complex(3,:,:)),':')
%plot(t, squeeze(Snfit_mag(3,:,:)),'--'), hold off
title('Alanine signals and fits')

return
%% Test fitting - 3-site

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(1:3,:,Iflips), TR, flips(1:3,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:2,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(1:3,:,Iflips)), TR, flips(1:3,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])
disp(['KPB  = ' num2str(kPB)])

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
k_fit = k_fit([[1:2],[7:8],[13:14]]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1:2],[7:8],[13:14]]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1:2],[7:8],[13:14]]); 
disp([['KPL  = ';'KPB  = '],num2str(reshape(k_fit,[2 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
subplot(223) , plot(t, squeeze(Sn(3,:,:)))
hold on, plot(t, squeeze(Snfit_complex(2,:,:)),':')
plot(t, squeeze(Snfit_mag(2,:,:)),'--'), hold off
title('Bicarb signals and fits')

%% Test fitting - 2-site

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
clear Sfit Snfit_complex Snfit_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est; params_fixed.R1B = R1B_est; params_fixed.R1A = R1A_est;
params_est.kPL = kPL_est; params_est.kPB = kPB_est; params_est.kPA = kPA_est;

for Iflips = 1:N_flip_schemes
    % no noise
    [params_fit(:,Iflips) Sfit(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Mxy(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(Sn(1:2,:,Iflips), TR, flips(1:2,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1,1:size(Mxy,2),  Iflips)] = fit_pyr_kinetics(abs(Sn(1:2,:,Iflips)), TR, flips(1:2,:,Iflips),params_fixed, params_est, std_noise, plot_fits);
end

disp('Input:')
disp(['KPL  = ' num2str(kPL)])

disp('Noiseless fit results:')
k_fit = struct2array(params_fit);
k_fit = k_fit([[1],[5],[9]]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])
disp('Noisy complex fit results:')
k_fit = struct2array(params_fitn_complex);
k_fit = k_fit([[1],[5],[9]]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])
disp('Noisy magnitude fit results:')
k_fit = struct2array(params_fitn_mag);
k_fit = k_fit([[1],[5],[9]]); 
disp([['KPL  = '],num2str(reshape(k_fit,[1 3]),2)])

figure
subplot(221) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(222) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(1,:,:)),':')
plot(t, squeeze(Snfit_mag(1,:,:)),'--'), hold off
title('Lactate signals and fits (dots=complex fit, dashed=magnitude)')
legend('constant','multiband', 'multiband variable flip')
