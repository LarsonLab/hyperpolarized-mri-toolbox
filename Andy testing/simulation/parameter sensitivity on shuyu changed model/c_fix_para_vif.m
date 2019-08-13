% ------------------------------------------------
% % simulation for vif
% ------------------------------------------------
clear all
disp('vif simulation');
% close all
% clc
tic

%--------------------
% flags
%--------------------
figure_flag=0;
pause_time=0;
plot_fits = 0;
%--------------------

% choose functions to test
sml_function = @changed_simulate_Nsite_perfusion_model;
fit_function = @changed_fit_kPL_perfusion;

% vif_array=[1,2,3,4,5,6,7,8,9,10]*0.1;
vif_array=linspace(0,1,20);
kpl_array=linspace(0.1,1,20);
% kpl_array=[2,4,6,8]*0.1;

% Test values
Tin = 0; Tacq = 20*3; TR = 3; N = Tacq/TR;
R1P = 1/25; R1L = 1/25; 
% kPL = 0.2;
% std_noise = 0.001;
kve = 0.02; vb = 0.09;

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

% Test over multiple combinations of flip angle schemes
flips(1:2,1:N,1) = repmat([10;30]*pi/180,[1 N]);
% flips(1,1:N,1) = [4,5,6,7,8,9,10,12,14,16,19,21,23,25,28,31,34,37,41,90]*pi/180;
flip_descripton{1} = 'constant, multi-band';
N_flip_schemes = size(flips,3);
t = [0:N-1]*TR + Tin;

% ------------------------------
% Choose to add noise or no
% ------------------------------
% noise_S = randn([2 N])*std_noise;  % same noise for all flip schedules
% noise_VIF = randn([2 N])*std_noise;  % same noise for all flip schedules
noise_S = zeros(2,N);  % same noise for all flip schedules
noise_VIF = 0;  % same noise for all flip schedules
% ------------------------------


for kpl_index=1:length(kpl_array) % kpl simulation
    kPL=kpl_array(kpl_index);

    for j = 1:length(vif_array) % vif simulation
        VIFscale = vif_array(j);
        % VIFscale = 0.5;
        % [Mxy_ev Mz_ev Mxy_iv Mz_iv] = sml_function...
        %     (Tin, [R1P R1L], [kPL 0], flips, TR, VIF(1,:), VIFscale, kve);

        [Mxy Mxy_iv Mxy_ev Mz Mz_iv Mz_ev] = sml_function...
        ([0;0], [R1P R1L], [kPL 0], flips, TR, VIF_z(1,:), VIFscale, kve);
        S(:,:) = Mxy;
        %----------------------------
        % cut the signal to be the same as the real data.
        %----------------------------
        % flips=flips(:,4:end,1);
        % t=t(4:end);
        % VIF=VIF(:,4:end);
        % Mxy_ev=Mxy_ev(:,4:end);
        % Mz_ev=Mz_ev(:,4:end);
        % Mxy_iv=Mxy_iv(:,4:end);
        % noise_S=noise_S(:,4:end);
        % Mz_iv=Mz_iv(:,4:end);
        %----------------------------

        %------------------------
        % Argue, also look into fit_kPL_perfused_voxel line 393, 394
        %------------------------
        % S(:,:) = Mxy_ev + Mxy_iv;
        % S(:,:) = Mxy_ev*(1-vb) + Mxy_iv*vb;
        %------------------------

        % Argue
        VIF_xy = Mxy_iv/VIFscale;
        % add noise
        Sn(:, : ) = S(:,:); % + noise_S;

        % initial parameter guesses
        R1P_est = R1P; R1L_est = R1L; kPL_est = kPL; 
         % kve_est = 0.02; vb_est = 0.1;
        kve_est = kve; 
        vb_est = vb;

        VIFscale_est= VIFscale;
        % VIFscale_est=VIFscale*0.1*(i-1);
        % VIFscale_est=VIFscale*10^(0.1*(i-11));
        % VIFscale_est=0.1*i;

        clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
        params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
        params_fixed.kve = kve_est; params_fixed.vb = vb_est;

        params_est.VIFscale = VIFscale_est;
        params_est.kPL = kPL_est;

        [params_fit(:) Sfit] = fit_function(S(:,:), VIF_xy(1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
        kpl_fit(kpl_index,j)=getfield(struct2table(params_fit),'kPL');
    end % vif simulation
end % kpl simulation
toc
%%
figure(2)
subplot(1,3,3)
hold on;
Legend=cell(length(kpl_array),1);
for kpl_index=1:length(kpl_array)
semilogx(vif_array,kpl_fit(kpl_index,:),'-o');
Legend{kpl_index}=['k_{pl} estimation = ',num2str(kpl_array(kpl_index))];
end
legend(Legend,'Location','northeastoutside');
title(['kPL (',num2str(kPL),') - VIF, fix all except kpl']);
xlabel('VIF');
ylabel('k_{PL}');
axis([0,1,0,1])

%%
figure
hold on;
plot(t,Sfit,'-*')
plot(t,Mxy_ev(2,:),'-o')
legend('fitted lactate','simulation lactate')