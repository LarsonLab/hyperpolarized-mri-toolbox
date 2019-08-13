% ------------------------------------------------
% % simulation
% ------------------------------------------------
clear all
disp('vif simulation');
% close all
% clc
tic
figure_flag=0;
% cut_position=2;
pause_time=0;
% choose fitting function to test

plot_fits = 0;

vif_array=[2,4,6,8,9,10]*0.1;
for j = 1:length(vif_array)
        % % Test values
        % Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
        % R1P = 1/25; R1L = 1/25; kPL = 0.5; std_noise = 0.001;
        % kve = kve_array(j); vb = 0.09;


    % Test values
    Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
    R1P = 1/25; R1L = 1/25; kPL = 0.5; std_noise = 0.001;
    kve = 0.02; vb = 0.09;
    VIFscale = vif_array(j);

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


    % Test over multiple combinations of flip angle schemes
    flips(1:2,1:N,1) = repmat([90;90]*pi/180,[1 N]);
    flips(1,1:N,1) = [4,5,6,7,8,9,10,12,14,16,19,21,23,25,28,31,34,37,41,90]*pi/180;
    flip_descripton{1} = 'constant, multi-band';
    N_flip_schemes = size(flips,3);
    t = [0:N-1]*TR + Tin;


    flip_description_array = [repmat('    ',N_flip_schemes,1),  char(flip_descripton(:))];

    % generate simulated data
    % noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
    % noise_VIF = randn([2 N])*std_noise;  % same noise for all flip schedules
    noise_S = zeros(2,N);  % same noise for all flip schedules
    noise_VIF = 0;  % same noise for all flip schedules

    [Mxy_ev Mz_ev Mxy_iv Mz_iv] = simulate_Nsite_perfused_voxel_model...
        (Mz0, [R1P R1L], [kPL 0], kve, vb, flips(:,:), TR, VIF*VIFscale);
    % cut the signal
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
    % VIF_xy_noisy = Mxy_iv + noise_VIF;
    VIF_xy = Mxy_iv/VIFscale;
    % add noise
    Sn(:, : ) = S(:,:);
     % + noise_S;

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
         % kve_est = 0.02; vb_est = 0.1;
        kve_est = kve; 
        vb_est = vb;
        % VIFscale_est=VIFscale*0.1*(i-1);
        VIFscale_est=VIFscale*10^(0.1*(i-11));
        % VIFscale_est=0.1*i;
        VIFscale_array(i)=VIFscale_est;
        clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
        params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
        
        params_fixed.kve = kve_est; params_fixed.vb = vb_est;


        params_est.VIFscale = VIFscale_est;
        params_est.kPL = kPL_est;
        [params_fit(:) Sfit(1:2,1:size(S,2))] = fit_function(S(:,:), VIF_xy, TR, flips(:,:), params_fixed, params_est, [], plot_fits);
        kpl_fit(j,i)=getfield(struct2table(params_fit),'kPL');
    end
end
toc
%%
kpl_vif_mean=mean(kpl_fit(1:14))
kpl_vif_std=std(kpl_fit(1:14))
% save('vif.mat','kpl_vif_mean','kpl_vif_std')
%%


figure(1)
subplot(1,3,3)
semilogx(VIFscale_array,kPL*ones(1,21));
hold on;
semilogx(VIFscale_array,kpl_fit,'-o');
% legend('True value of kPL','k_{PL} estimated')
legend('True value of kPL','vif =0.2','vif=0.4','vif=0.6','vif=0.8','vif=0.9','vif=1')

title(['kPL (',num2str(kPL),') - VIF, fix others']);
xlabel('VIF error factor');
ylabel('k_{PL}');
axis([0,10,0.46,0.58])
