 % This file is to generate the plot for 
% kpl - noise standard deviation.

% ------------------------------------------------
% % simulation
% ------------------------------------------------
clear all
close all
clc
tic
figure_flag=0;
cut_position=2;
pause_time=0;
% choose fitting function to test

plot_fits = 0;

% Test values
Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; kPL = 0.2; 
% std_noise = 0.001;
kve = 0.05; vb = 0.9;

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
VIF_z = [input_function; zeros(1,N)];
VIFscale = 0.7;

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = repmat([20;30]*pi/180,[1 N]);
% flips(1,1:N,1) = [4,5,6,7,8,9,10,12,14,16,19,21,23,25,28,31,34,37,41,90]*pi/180;
flip_descripton{1} = 'constant, multi-band';
N_flip_schemes = size(flips,3);
t = [0:N-1]*TR + Tin;

if figure_flag==1
figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi)
title('Pyruvate flips')
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi)
title('Lactate flips')
legend(flip_descripton)
drawnow, pause(pause_time)
end

flip_description_array = [repmat('    ',N_flip_schemes,1),  char(flip_descripton(:))];

% generate simulated data
std_noise_base = 0.001;
M=100; % repeat noise
R=10; % increase std times
for j=1:R
    j
    std_noise=std_noise_base*j;
    std_noise_array(j)=std_noise;
    noise_S = randn([2*M N])*std_noise;
    noise_VIF = randn([2*M N])*std_noise;

    % [Mxy_ev Mz_ev Mxy_iv Mz_iv] = simulate_Nsite_perfused_voxel_model...
    %     (Mz0, [R1P R1L], [kPL 0], kve, vb, flips(:,:), TR, VIF*VIFscale);
[Mxy Mxy_iv Mxy_ev Mz Mz_iv Mz_ev] = changed_simulate_Nsite_perfusion_model...
        ([0;0], [R1P R1L], [kPL 0], flips, TR, VIF_z(1,:), VIFscale, kve);
    S(:,:) = Mxy;
    VIF_xy_noisy = Mxy_iv + noise_VIF;
    % VIF_xy = Mxy_iv/VIFscale;
    VIF_xy = [Mxy_iv/VIFscale;zeros(1,20)];
    % add noise
    % Sn(:, : ) = S(:,:) + noise_S;



    % S(:,:) = Mxy_ev + Mxy_iv;
    % VIF_xy = Mxy_iv;
    % VIF_xy_noisy = Mxy_iv + noise_VIF;
    % add noise
    S_rep=repmat(S,M,1);
    VIF_xy_rep=repmat(VIF_xy,M,1);
    VIF_xy_noisy=VIF_xy_rep+noise_VIF;
    Sn(:, : ) = S_rep(:,:) + noise_S;


    % clearvars -except t S_rep Sn S
    Sn_max=[max(Sn(1:2:end,:));max(Sn(2:2:end,:))];
    Sn_min=[min(Sn(1:2:end,:));min(Sn(2:2:end,:))];
    disp(['estimate vif from simulation: ',num2str(sum(S(1,1))/sum(VIF_xy(1,1)))]);
    VIFscale_est=sum(S(1,1))/sum(VIF_xy(1,1))
    % figure
    % hold on;
    % plot(t,S);
    % plot(t,Sn_min,'--',t,Sn_max,'--');
    % legend('no noise pyr','no noise lac','noisy pyr','noisy lac');
    % title('simulation result');


    disp('------------------------');
    disp('new model');
    fit_function = @changed_fit_kPL_perfusion;

    % initial parameter guesses


    R1P_est = 1/25; R1L_est = 1/25; kPL_est = .2; 
     % kve_est = 0.02; vb_est = 0.1;
    kve_est = kve; 
    vb_est = vb;
    % VIFscale_est=VIFscale;

    clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
    params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;

    params_fixed.kve = kve_est;
     % params_fixed.vb = vb_est;
    params_fixed.VIFscale = VIFscale_est;
    params_est.kPL = kPL_est;
    % add IV function
    % [params_fit_optimal(:) Sfit(1:2,1:size(S,2))] = fit_function(S(:,:), VIF_xy, TR, flips(:,:), params_fixed, params_est, [], plot_fits);
    [params_fit_optimal(:) Sfit2] = fit_function(S(:,:), VIF_xy(1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
    for i=1:M
    % [params_fit_noisy(i) Sfit_noisy(2*i-1:2*i,1:size(S,2))] = fit_function(Sn(2*i-1:2*i,:), VIF_xy_rep(2*i-1:2*i,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
    [params_fit_noisy(i) Sfit2(i,1:size(S,2))] = fit_function(Sn(2*i-1:2*i,:), VIF_xy_rep(2*i-1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
    end
    %%
    % disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, kPL))
    % disp(sprintf('Input :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL, VIFscale, kve,vb))
    % disp(sprintf('Estimate :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL_est, VIFscale_est, kve_est,vb_est))
    % disp('Noiseless fit results:')
    % disp(['kPL  = ',num2str(getfield(struct2table(params_fit_optimal),'kPL'))]);fprintf('%c',  8);
    % disp(['  ','VIFscale  = ',num2str(getfield(struct2table(params_fit_optimal),'VIFscale'))]);fprintf('%c',  8);
    % disp(['  ','kve  = ',num2str(getfield(struct2table(params_fit_optimal),'kve'))]);fprintf('%c',  8);
    % disp(['  ','vb  = ',num2str(getfield(struct2table(params_fit_optimal),'vb'))]);
    % disp('VIFscale  = ']); disp([num2str(getfield(struct2table(params_fit),'VIFscale'))])
    kpl_optimal(j)=getfield(struct2table(params_fit_optimal),'kPL');
    kpl_std(j)=std(getfield(struct2table(params_fit_noisy),'kPL'));
    kpl_mean(j)=mean(getfield(struct2table(params_fit_noisy),'kPL'));
    maxlac(j)=max(Sfit2(:));
    % kpl_min(j)=min(getfield(struct2table(params_fit_noisy),'kPL'));
end
toc

%%
% kpl_mean2=kpl_mean;
% kpl_std2=kpl_std;
% 
figure
hold on;
% subplot(1,3,1)
plot(maxlac./std_noise_array,kpl_optimal);
errorbar(maxlac./std_noise_array,kpl_mean,kpl_std,'-o');
legend('True value of kPL','Mean and std')
xlabel('std for noise (both signal and VIF)');
ylabel('k_{PL}');
title('kPL - noise (flip: 10,30)');
% axis([1e-3,1.1*10e-3,0,.5])
% subplot(1,3,2)
% plot(std_noise_array,kpl_optimal);
% errorbar(std_noise_array,kpl_mean2,kpl_std2,'-o');
% title('kPL - noise (flip: 90,90)');
% legend('True value of kPL','Mean and std')
% axis([1e-3,1.1*10e-3,-1.5,3.5])
% xlabel('std for noise (both signal and VIF)');
% ylabel('k_{PL}');
% subplot(1,3,3)
% plot(std_noise_array,kpl_optimal);
% errorbar(std_noise_array,kpl_mean1,kpl_std1,'-*');
% legend('True value of kPL','Mean and std')

% axis([1e-3,1.1*10e-3,-1.5,3.5])
% % legend('True value of kPL','Mean and std - 10 flip','Mean and std - 90 flip','mean and std - data flip')
% title('kPL - noise (flip: data flip)');
% xlabel('std for noise (both signal and VIF)');
% ylabel('k_{PL}');