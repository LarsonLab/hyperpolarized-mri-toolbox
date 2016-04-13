% Script for testing fit_kinetics_* kinetic model fitting functions

clear all

% Test values
Tin = 6; Tacq = 30; TR = 3; N = Tacq/TR;
R1 = [1/25 1/25]; KPL = 0.05; std_noise = 0.01;
k12 = 0.05; Minit = [1; 0.07];

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = ones(2,N)*30*pi/180;  % constant, single-band
flips(1:2,1:N,2) = repmat([20;35]*pi/180,[1 N]);  % constant, multi-band
flips(1:2,1:N,3) = repmat(vfa_const_amp(N, pi/2), [2,1]);  % RF compensated variable flip angle
flips(1:2,1:N,4) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % T1-effective variable flip angle
		    vfa_const_amp(N, pi/2, exp(-TR * ( - k12)))];
flips(1:2,1:N,5) = [vfa_const_amp(N, pi/2, exp(-TR * ( k12))); ... % max lactate SNR variable flip angle
		    vfa_opt_signal(N, exp(-TR * ( R1(2))))];
flips(1:2,1:N,6) = [vfa_const_amp(N, pi/2, exp(-TR * (k12))); ... % saturation recovery
		    ones(1,N)*pi/2];

         t = [0:N-1]*TR + Tin;

        figure
 subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
 subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Lactate flips')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')
        
%% Test fitting
 for Iflips = 1:size(flips,3)
     [Mxy Mz] = simulate_2site_model(Tin, R1, [KPL 0], flips(:,:,Iflips), TR);
     % no noise
    [R1fit(Iflips) KPLfit(Iflips) Sfit(1:size(Mxy,1), 1:size(Mxy,2),  Iflips)] = fit_kinetics_vfa_sameT1(Mxy, TR, flips(:,:,Iflips));
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2),  Iflips) = Mxy + randn(size(Mxy))*std_noise;
    [R1fitn(Iflips) KPLfitn(Iflips) Snfit(1:size(Mxy,1), 1:size(Mxy,2),  Iflips)] = fit_kinetics_vfa_sameT1(Sn(:,:,Iflips), TR, flips(:,:,Iflips));
 end

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1(1), R1(2), KPL))
disp('Noiseless fit results:')
disp(['R1fit  = ' num2str(R1fit)])
disp(['KPLfit = ' num2str(KPLfit)])
disp('Noisy fit results:')
disp(['R1fitn  = ' num2str(R1fitn)])
disp(['KPLfitn = ' num2str(KPLfitn)])
 
  figure
 subplot(121) , plot(t, squeeze(Sn(1,:,:)))
    hold on, plot(t, squeeze(Snfit(1,:,:)),'--')
title('Pyruvate signals and fits (dashed)')
 subplot(122) , plot(t, squeeze(Sn(2,:,:)))
    hold on, plot(t, squeeze(Snfit(2,:,:)),'--')
title('Lactate signals and fits (dashed)')
legend('constant','multiband', 'vfa', 'T1-effective vfa', 'max lactate SNR vfa', 'Saturation Recovery')
