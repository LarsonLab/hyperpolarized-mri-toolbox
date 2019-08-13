% SNR=peak(lac)/std_nosie=5,10,20,50,100

% ------------------------------------------------
% % simulation
% ------------------------------------------------
clear all
close all
clc
figure_flag=1;
cut_position=2;
pause_time=0;
% choose fitting function to test

plot_fits = 0;

% Test values
Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; kPL = 0.2; std_noise = 0.001;
kve = 0.02; vb = 0.2;

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
VIFscale = 0.9;

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = repmat([20;30]*pi/180,[1 N]);
% flips(1,1:N,1) = [4,5,6,7,8,9,10,12,14,16,19,21,23,25,28,31,34,37,41,90]*pi/180;

flip_descripton{1} = 'constant, multi-band';
N_flip_schemes = size(flips,3);
t = [0:N-1]*TR + Tin;
if figure_flag==1
figure
subplot(121) , plot(t, squeeze(flips(1,:,:))*180/pi,'-o')
title('Pyruvate flips')
axis([0,60,0,90])
subplot(122) , plot(t, squeeze(flips(2,:,:))*180/pi,'-o')
title('Lactate flips')
legend(flip_descripton)
axis([0,60,0,90])
drawnow, pause(pause_time)
end

flip_description_array = [repmat('    ',N_flip_schemes,1),  char(flip_descripton(:))];

% generate simulated data
% noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
% noise_VIF = randn([2 N])*std_noise;  % same noise for all flip schedules
noise_S = zeros(2,N);  % same noise for all flip schedules
noise_VIF = 0;  % same noise for all flip schedules
    [Mxy Mxy_iv Mxy_ev Mz Mz_iv Mz_ev] = changed_simulate_Nsite_perfusion_model...
        ([0;0], [R1P R1L], [kPL 0], flips, TR, VIF_z(1,:), VIFscale, kve);
    S(:,:) = Mxy;
    VIF_xy_noisy = Mxy_iv + noise_VIF;
    % VIF_xy = Mxy_iv/VIFscale;
    VIF_xy = Mxy_iv/VIFscale;
    % add noise
    Sn(:, : ) = S(:,:) + noise_S;
if figure_flag==1
    figure
    subplot(221)
    plot(t, Mz_iv), title('M_Z intravascular')
    legend('pyruvate','lactate')
    subplot(222)
    plot(t, Mz_ev), title('M_Z extravascular')
    legend('pyruvate','lactate')
    subplot(223)
    plot(t, Mxy_iv), title('M_{XY} intravascular')
    legend('pyruvate','lactate')
    subplot(224)
    plot(t, Mxy_ev), title('M_{XY} extravascular')
    legend('pyruvate','lactate')
    drawnow, pause(pause_time)
end
figure(199)
hold on;
plot(S(1,:),'-s')
plot(VIF_xy(1,:),'-o')
legend('s','vif xy')
hold off;
%%
disp(['estimate vif from simulation: ',num2str(sum(S(1,1))/sum(VIF_xy(1,1)))]);
% ------------------------------------------------
% % old model
% ------------------------------------------------
disp('------------------------');
disp('old model');

fit_function = @fit_kPL;
% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; kPL_est = .02;

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;
% params_est.L0_start=1;

[params_fit(:) Sfit1(1:size(S,2))] = fit_function(S(:,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);

disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, kPL))
disp('Noiseless fit results:')
disp(['kPL  = ']); disp([num2str(struct2array(params_fit).'), flip_description_array])

titles={'simulation','fitting'};
if figure_flag==1
figure(99)
subplot(121) ,hold on; plot(t, squeeze(S(1,:)),'-s')
% hold on, plot(t, squeeze(Sfit1(1,:)),'s')
% legend(char(titles))
title('Pyruvate signals')
subplot(122) , plot(t, squeeze(S(2,:)),'-*')
hold on, plot(t, squeeze(Sfit1(1,:)),'--s')
% legend(char(titles))
% title('old Lactate signals')
% drawnow, pause(0.5)
end

% ------------------------------------------------
% % perfused model
% ------------------------------------------------
disp('------------------------');
disp('new model');
fit_function = @changed_fit_kPL_perfusion;

% initial parameter guesses
R1P_est = 1/25; R1L_est = 1/25; kPL_est = .2; 
% kve_est = 0.02; vb_est = 0.1;
kve_est = kve; 
vb_est = vb;
VIFscale_est=VIFscale;

%% Test fitting - fit kPL only
% disp('Fitting kPL, with fixed relaxation rates:')
% disp('Fixing relaxation rates improves the precision of fitting, but potential')
% disp('for bias in fits when incorrect relaxation rate is used')
% disp(' ')

clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
% params_est={};
% params_fixed.kve = kve_est; params_fixed.vb = vb_est;
params_est.VIFscale = VIFscale_est;
params_est.kPL = kPL_est;
% add IV function
[params_fit(:) Sfit2] = fit_function(S(:,:), VIF_xy(1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
%%
disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, kPL))
disp(sprintf('Input :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL, VIFscale, kve,vb))
disp(sprintf('Estimate :\n kPL = %f VIFscale= %f kve= %f vb=%f', kPL_est, VIFscale_est, kve_est,vb_est))
disp('Noiseless fit results:')
disp(['kPL  = ',num2str(getfield(struct2table(params_fit),'kPL'))]);fprintf('%c',  8);
disp(['  ','VIFscale  = ',num2str(getfield(struct2table(params_fit),'VIFScale'))]);fprintf('%c',  8);
% disp(['  ','kve  = ',num2str(getfield(struct2table(params_fit),'kve'))]);fprintf('%c',  8);
% disp(['  ','vb  = ',num2str(getfield(struct2table(params_fit),'vb'))]);
% disp('VIFscale  = ']); disp([num2str(getfield(struct2table(params_fit),'VIFscale'))])






if figure_flag==1
figure(99)
subplot(121)
hold on;
% plot(t, squeeze(Sfit2(1,:)),'--o')
% plot(t,fitvif,':*')
% plot(t,Mxy_ev(1,:),':')
% plot(t,Mxy_iv(1,:),':')
% legend('simulation-truth','old model fitted','new model fitted','fitted VIF')
legend('simulation-truth')
% ,'ev','iv')
title('pyruvate comparison')
subplot(122)
hold on;
plot(t, squeeze(Sfit2(1,:)),'--o')

% plot(t,Mxy_ev(2,:),':')
% plot(t,Mxy_iv(2,:),':')
legend('simulation truth','old fitted','new fitted')
% ,'ev','iv')
title('lactate comparison')
drawnow, pause(pause_time)

figure
subplot(121)
hold on;
plot(t, squeeze(S(1,:)),'--o')

plot(t,Mxy_ev(1,:),'--*')
plot(t,Mxy_iv(1,:),'--s')
legend('simulation-truth','ev','iv')
title('pyruvate comparison')
subplot(122)
hold on;
plot(t, squeeze(S(2,:)),'--o')

plot(t,Mxy_ev(2,:),'--*')
plot(t,Mxy_iv(1,:),'--s')
legend('simulation truth','ev','iv')
title('lactate comparison')
drawnow, pause(pause_time)

% titles={'simulation','fitting'};
% figure
% subplot(121) , plot(t, squeeze(S(1,:)),'-*')
% hold on, plot(t, squeeze(Sfit2(1,:)),'s')
% legend(char(titles))
% title('Pyruvate signals')
% subplot(122) , plot(t, squeeze(S(2,:)),'-*')
% hold on, plot(t, squeeze(Sfit2(2,:)),'s')
% legend(char(titles))
% title('new Lactate signals')
% drawnow, pause(pause_time)

end
%%
figure
hold on;
plot(VIF_z(1,:))
plot(VIF_xy(1,:))
legend('vif z','vif xy')
hold off;