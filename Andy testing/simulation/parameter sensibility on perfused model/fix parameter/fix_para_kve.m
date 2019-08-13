% ------------------------------------------------
% % simulation
% ------------------------------------------------
clear all
disp('kve simulation');
% close all
% clc
tic
figure_flag=0;
% cut_position=2;
pause_time=0;
% choose fitting function to test

plot_fits = 0;
kve_array=[0.005,0.01,0.02,0.03,0.04];
for j = 1:length(kve_array)
    % Test values
    Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
    R1P = 1/25; R1L = 1/25; kPL = 0.5; std_noise = 0.001;
    kve = kve_array(j); vb = 0.09;

    input_function = zeros(1,N);

    input_condition = 1; % choose from various simulated starting conditions
    switch input_condition
        case 1
            % gamma variate input function - most realistic
            t = [1:N]*TR;
            Tarrival = 0;  Tbolus = 12;
            input_function = realistic_input_function(N, TR, Tarrival, Tbolus);
            Mz0 = [0,0]; 
        case 2
            % boxcar input function
            Tbolus = 12;  Tarrival = 0;
            Ibolus = [1:round(Tbolus/TR)] + round(Tarrival/TR);
            Rinj = 1/Tbolus; % normalize for a total magnetization input = 1
            Mz0 = [0,0]; input_function(Ibolus) =  Rinj*TR;

    end
    VIF = [input_function; zeros(1,N)];
    VIFscale = 1;

    % Test over multiple combinations of flip angle schemes
    flips(1:2,1:N,1) = repmat([90;90]*pi/180,[1 N]);
    flips(1,1:N,1) = [4,5,6,7,8,9,10,12,14,16,19,21,23,25,28,31,34,37,41,90]*pi/180;
    flip_descripton{1} = 'constant, multi-band';
    N_flip_schemes = size(flips,3);
    t = [0:N-1]*TR + Tin;

    noise_S = zeros(2,N);  % same noise for all flip schedules
    noise_VIF = 0;  % same noise for all flip schedules

    [Mxy_ev Mz_ev Mxy_iv Mz_iv] = simulate_Nsite_perfused_voxel_model...
        (Mz0, [R1P R1L], [kPL 0], kve, vb, flips(:,:), TR, VIF*VIFscale);
    flips=flips(:,4:end,1);
    t=t(4:end);
    VIF=VIF(:,4:end);
    Mxy_ev=Mxy_ev(:,4:end);
    Mz_ev=Mz_ev(:,4:end);
    Mxy_iv=Mxy_iv(:,4:end);
    noise_S=noise_S(:,4:end);
    Mz_iv=Mz_iv(:,4:end);
    % Argue
    S(:,:) = Mxy_ev + Mxy_iv;
    % S(:,:) = Mxy_ev*(1-vb) + Mxy_iv*vb;

    % Argue
    % VIF xy construction
    VIF_xy_noisy = Mxy_iv + noise_VIF;
    % VIF_xy = Mxy_iv;
    VIF_xy = Mxy_iv/VIFscale;
    % add noise
    Sn(:, : ) = S(:,:) + noise_S;

    % ------------------------------------------------
    % % perfused model
    % ------------------------------------------------
    % disp('------------------------');
    % disp('new model');
    fit_function = @fit_kPL_perfused_voxel;

    % initial parameter guesses
    % Test values
    % Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
    % R1P = 1/25; R1L = 1/25; kPL = 0.2; std_noise = 0.001;
    % kve = 0.05; vb = 0.9;

    for i = 1:21

        R1P_est = R1P; R1L_est = R1L; kPL_est = kPL; 
        vb_est = vb;
        VIFscale_est=VIFscale;
        
        kve_est = kve*0.1*(i-1); 
        kve_factor_array(i)=0.1*(i-1);
        clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
        params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
        params_fixed.vb = vb_est;
        params_fixed.VIFscale = VIFscale_est;


        params_est.kPL = kPL_est;
        params_fixed.kve = kve_est;
        [params_fit(:) Sfit(1:2,1:size(S,2))] = fit_function(S(:,:), VIF_xy, TR, flips(:,:), params_fixed, params_est, [], plot_fits);
        kpl_fit(j,i)=getfield(struct2table(params_fit),'kPL');
%         disp(sprintf('Input :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL, VIFscale, kve,vb))
% disp(sprintf('Estimate :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL_est, VIFscale_est, kve_est,vb_est))
% disp('Noiseless fit results:')
% disp(['kPL  = ',num2str(getfield(struct2table(params_fit),'kPL'))]);fprintf('%c',  8);
% disp(['  ','VIFscale  = ',num2str(getfield(struct2table(params_fit),'VIFscale'))]);fprintf('%c',  8);
% disp(['  ','kve  = ',num2str(getfield(struct2table(params_fit),'kve'))]);fprintf('%c',  8);
% disp(['  ','vb  = ',num2str(getfield(struct2table(params_fit),'vb'))]);
% disp('-------------');
    end
end % kve simulation
toc

kpl_kve_mean=mean(kpl_fit);
kpl_kve_std=std(kpl_fit);
% save('kve.mat','kpl_kve_mean','kpl_kve_std')
%%


figure(1)
subplot(1,3,2)
plot(kve_factor_array,kPL*ones(1,21));
hold on;
plot(kve_factor_array,kpl_fit,'--o');
legend('True value of kPL','k_{ve} =0.005','k_{ve}=0.01','k_{ve}=0.02','k_{ve}=0.03','k_{ve}=0.04')
 % title('kPL(0.1) - k_{ve}');
axis([0,2,0.1,0.8])
title(['kPL (',num2str(kPL),') - k_{ve}, fix all except kpl']);
xlabel('k_{ve} error factor');
ylabel('k_{PL}');

