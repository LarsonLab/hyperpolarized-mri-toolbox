clear all
NMC = 100;  % less for quicker testing

% default experiment values
experiment.R1P = 1/25;  experiment.R1L =1/25;  experiment.kPL = 0.02; experiment.std_noise = 0.01;
experiment.Tarrival = 4;  experiment.Tbolus = 8;

for  est_R1L = 0
    for fit_model = 1:3
        disp('Running Monte Carlo Simulation')
        fit_description = [];
        
        clear params_est params_fixed acq fitting
        
        % default fitting parameters
        R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;
        params_fixed.R1P = R1P_est;
        params_est.kPL = kPL_est;
        if est_R1L
            fit_description = [fit_description, '  Fitting T1L'];
            params_est.R1L = R1L_est;
        else
            fit_description = [fit_description, '  Fixed T1L'];
            params_fixed.R1L = R1L_est;
        end
        
%        params_fixed.L0_start = 0;  % ok to be free parameter?

        % allowing for fit of R1L increases variability substantially
        % R1P minimal change
        % constraining R1L to narrow range maybe reasonable compromise
        
        switch fit_model
            case 1
           fit_description = [fit_description, '  Inputless fitting (forward)'];
           fitting.fit_fcn = @fit_kPL;
            case 2
           fit_description = [fit_description, '  Inputless fitting (reverse)'];
           fitting.fit_fcn = @fit_kPL_reverse;
            case 3
           fit_description = [fit_description, '  Inputless fitting (forward+reverse)'];
           fitting.fit_fcn = @fit_kPL_bidirectional;
            case 4 % start from max?
        end
        
        
        
        % 2D dynamic 10/20 flips
        Tacq = 90; acq.TR = 5; acq.N = Tacq/acq.TR;
        Npe = 8; Nall = acq.N * Npe;
       acq.flips(1:2,1:acq.N) = repmat(acos(cos([10*pi/180;
       20*pi/180]).^Npe), [1 acq.N]);  % reverse failing sometimes

% Tacq = 45; acq.TR = 3; acq.N = Tacq/acq.TR;
%         acq.flips = [vfa_const_amp(acq.N, pi/2, exp(-acq.TR * ( 0.05))); ... % max lactate SNR variable flip angle
%     vfa_opt_signal(acq.N, exp(-acq.TR * ( experiment.R1L)))];
% 
         %                   acq.flips = repmat([10*pi/180; 40*pi/180], [1 acq.N]);  % reverse failing sometimes
        
        
        fitting.params_est = params_est; fitting.params_fixed = params_fixed;
        fitting.NMC = NMC;
        
        disp(fit_description)
        [results, hdata, hsim ] = HP_montecarlo_evaluation( acq, fitting, experiment );
        hdata.Name = fit_description; hsim.Name =fit_description;
        
    end
end

return

%%

N = acq.N; TR = acq.TR; std_noise = experiment.std_noise;
input_function = realistic_input_function(acq.N, acq.TR, experiment.Tarrival, experiment.Tbolus);
Mz0 = [0,0]';

%         acq.flips(1:2,1:acq.N) = repmat(acos(cos([10*pi/180; 20*pi/180]).^Npe), [1 acq.N]);
flips = acq.flips;
% flips(1:2,1:N) = ones(2,N)*30*pi/180;  % constant, single-band
% flips(1:2,1:N) = repmat([20;35]*pi/180,[1 N]);  % constant, multi-band
% flips(1:2,1:N) = [vfa_const_amp(N, pi/2, exp(-TR * ( 0.05))); ... % max lactate SNR variable flip angle
%     vfa_opt_signal(N, exp(-TR * ( experiment.R1L)))];
% 
% reverse failing for high flip angles


% generate simulated data
noise_S = randn([2 N])*experiment.std_noise;  % same noise for all flip schedules

    [Mxy(1:2, 1:N), Mz] = simulate_2site_model(Mz0, [experiment.R1P experiment.R1L], [experiment.kPL 0], flips, TR, input_function);
%    Mxy = Mxy(:,1:2:end) + Mxy(:,2:2:end);
    % add noise
    Sn(1:size(Mxy,1), 1:size(Mxy,2)) = Mxy(:,:) + noise_S;

    clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

    
    plot_fits = 1;

    Iflips = 1;
    [params_fit(:,Iflips) Sfit(1:size(Mxy,2),  Iflips)] = fit_kPL_bidirectional(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);

    return
    % no noise
    [params_fit(:,Iflips) Sfit(1:size(Mxy,2),  Iflips)] = fit_kPL(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex(:,Iflips) Snfit_complex(1:size(Mxy,2),  Iflips)] = fit_kPL(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag(:,Iflips) Snfit_mag(1:size(Mxy,2),  Iflips)] = fit_kPL(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);

    % no noise
    [params_fit_reverse(:,Iflips) Sfit_reverse(1:size(Mxy,2),  Iflips)] = fit_kPL_reverse(Mxy(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % add noise
    [params_fitn_complex_reverse(:,Iflips) Snfit_complex_reverse(1:size(Mxy,2),  Iflips)] = fit_kPL_reverse(Sn(:,:,Iflips), TR, flips(:,:,Iflips), params_fixed, params_est, [], plot_fits);
    
    % magnitude fitting with noise
    [params_fitn_mag_reverse(:,Iflips) Snfit_mag_reverse(1:size(Mxy,2),  Iflips)] = fit_kPL_reverse(abs(Sn(:,:,Iflips)), TR, flips(:,:,Iflips),params_fixed, params_est, std_noise, plot_fits);

    
    disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', experiment.R1P, experiment.R1L, experiment.kPL))
disp('*Forward fits*')
disp('Noiseless fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fit).'))
disp('Noisy complex fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_complex).'))
disp('Noisy magnitude fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_mag).'))

disp('*Reverse fits*')
disp('Noiseless fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fit_reverse).'))
disp('Noisy complex fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_complex_reverse).'))
disp('Noisy magnitude fit results:')
disp(['KPL  = ']); disp(num2str(struct2array(params_fitn_mag_reverse).'))

t = [0:N-1]*TR;

figure
subplot(131) , plot(t, squeeze(Sn(1,:,:)))
title('Pyruvate signals')
subplot(132) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex(:,:)),':')
plot(t, squeeze(Snfit_mag(:,:)),'--')
title('kPL fit: Lactate signals and fits (dots=complex fit, dashed=magnitude)')
%legend('constant','multiband', 'multiband variable flip')
subplot(133) , plot(t, squeeze(Sn(2,:,:)))
hold on, plot(t, squeeze(Snfit_complex_reverse(:,:)),':')
plot(t, squeeze(Snfit_mag_reverse(:,:)),'--')
title('kPL fit reverse: Lactate signals and fits (dots=complex fit, dashed=magnitude)')

