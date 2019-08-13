% Andy's comment:
% 	This file is a modified version of Prof. Larson's test_fit_kPL_fcn.m
%   This file is for testing over the data collected by Renuka Sriram.
% 	Thanks for all the help from Shuyu Tang
% 	Script for testing fit_kPL kinetic model fitting functions

close all
clear all;
clc;

% data_dir = './Andy testing/Andy_data';
data_dir = './Andy testing/others/debug_data';
fitted_dir='./Andy testing/others/Andy_fitted';
tumors={'/786O','/A498','/uok262'};
tumors_num=length(tumors);
N = 20; % should double check with Renuka Sriram
fit_function = @fit_kPL;
line = 4; lac_sheet = 1; flip_sheet = 4;
plot_flag=3; plot_fits = 0; plot_mag = 0; load_avg=1; save_flag=1; disp_1xlsx_flag=0;

% picking data
thresh_above=3; % number of points above noise mean
noise_multiplyer=1.8; % multiplyer to the noise mean

for A_scaler=0.5:0.1:0.7 % loop for one tumor folder

	tumor_type='';
	[filenames,files_num] = A_get_filenames([data_dir,tumor_type]);
	filenames=char(filenames);
	% tumor_kpl=zeros(6,files_num);
	tumor_kpl=zeros(4,files_num);
	errors=zeros(1,files_num);
	% if mod(plot_flag,10)==3
	% 	figure
	% end

	% --------------------
	% load data;
	% --------------------
	for file_counter=1:files_num % loop for one xlsx file
		filename = filenames(file_counter,:);
		% filename='JW342_A498_122016_spspdyn copy.xlsx'
		filename=filename(1:strfind(filename, '.xlsx')-1);
		dir=[data_dir,tumor_type,'/',filename,'.xlsx'];
		if load_avg==1
			[flips,Mxy,t,TR,Tin,std_noise,mean_noise] = A_load_avg_data(dir,flip_sheet,N);
	    else
			[flips,Mxy,t,TR,std_noise] = A_load_vexel_data(dir,line,lac_sheet,flip_sheet,N);
	    end
	    % --------------------
		% Post process for input data
		% --------------------
	    std_noise=std_noise*sqrt(pi/2);
	    % first three flip is zero
	    Mxy=Mxy(:,4:end,:);
	    flips=flips(:,4:end);
	    t=t(4:end);
		% Test values;
		R1P = 1 / 25; R1L = 1 / 25; kPL = 0.05;
		input_function = zeros(1, N);

		mean_noise=ones(size(t))*mean_noise;
		use_flag=nnz(Mxy(1,:,1)>mean_noise*noise_multiplyer)>thresh_above;
		if use_flag==1
			use_flag=nnz(Mxy(1,:,2)>mean_noise*noise_multiplyer)>thresh_above;
		end
		% use_flag
		% --------------------
		% plot original data;
		% --------------------
		% plot flip angles;
		if plot_flag/10>=1
			figure
			subplot(221) , plot(t, squeeze(flips(1, : )) * 180 / pi)
			title('Pyruvate flips')
			subplot(222) , plot(t, squeeze(flips(2, : )) * 180 / pi,'o');
			hold on; plot(t, squeeze(flips(2, : )) * 180 / pi)
			title('Lactate flips')
			% plot pyr and lac
			subplot(2,2,3:4) , plot(t, squeeze(Mxy(1, : , 1)),'r')
			hold on
			plot(t, squeeze(Mxy(2, :,1 )),'r--')
			plot(t, squeeze(Mxy(1, :,2 )),'b')
			plot(t, squeeze(Mxy(2, :,2 )),'b--')
			plot(t, squeeze(Mxy(1, :,3 )),'g')
			plot(t, squeeze(Mxy(2, :,3 )),'g--')
			legend('Pyruvate tumor','Lactate tumor','Pyruvate control','Lactate control','Pyruvate noise','Lactate noise')
			title(['Mxy data',filename(1:10)])

		end

		if use_flag==1
			% initial parameter guesses;
			R1P_est = 1 / 25; R1L_est = 1 / 25; kPL_est = .02;
			% --------------------
			% Test fitting - fit kPL only;
			% --------------------
			if disp_1xlsx_flag==1
				disp('Fitting kPL, with fixed relaxation rates:')
			end
			
			% --------------------
			% Fitting with noise 1
			% --------------------

			clear params_fixed_1 params_est_1;
			clear params_fit1_t params_fit1_c params_fit1_n;
			clear sig_fit1_t sig_fit1_c sig_fit1_n;
			params_fixed_1.R1P = R1P_est; params_fixed_1.R1L = R1L_est;
			params_est_1.kPL = kPL_est;
			fig=figure
			subplot(2,3,1)
			% fitting for tumor;
			[params_fit1_t(:) sig_fit1_t(1:size(Mxy, 2))] = fit_function(Mxy(:, :,1), TR, flips( :, :), params_fixed_1, params_est_1, std_noise, plot_fits,A_scaler);
			title([num2str(A_scaler),', nml tmr, kpl: ',num2str(struct2array(params_fit1_t))])
			% fitting for control;
			
			subplot(2,3,2)
			% Mxy(2, :,2)=Mxy(2, :,2)*2;
			[params_fit1_c(:) sig_fit1_c(1:size(Mxy, 2))] = fit_function(Mxy(:, :,2), TR, flips( :, :), params_fixed_1, params_est_1, std_noise, plot_fits,A_scaler);
			title(['nml cntl, kpl: ',num2str(struct2array(params_fit1_c))])
			% fitting for noise;
			% [params_fit1_n(:) sig_fit1_n(1:size(Mxy, 2))] = fit_function(Mxy(:, :,3), TR, flips( :, :), params_fixed_1, params_est_1, std_noise, plot_fits);

			% --------------------
			% Fitting without noise 2
			% --------------------
			clear params_fixed_2 params_est_2;
			clear params_fit2_t params_fit2_c params_fit2_n;
			clear sig_fit2_t sig_fit2_c sig_fit2_n;
			params_fixed_2.R1P = R1P_est; params_fixed_2.R1L = R1L_est;
			params_est_2.kPL = kPL_est;
			% fitting for tumor;
			subplot(2,3,4)
			[params_fit2_t(:) sig_fit2_t(1:size(Mxy, 2))] = fit_function(Mxy(:, :,1), TR, flips( :, :), params_fixed_2, params_est_2, [], plot_fits,A_scaler);
			title([' ls tmr, kpl: ',num2str(struct2array(params_fit2_t))])
			% fitting for control;
			subplot(2,3,5)
			[params_fit2_c(:) sig_fit2_c(1:size(Mxy, 2))] = fit_function(Mxy(:, :,2), TR, flips( :, :), params_fixed_2, params_est_2, [], plot_fits,A_scaler);
			title([' ls cntl, kpl: ',num2str(struct2array(params_fit2_c))])
			% fitting for noise;
			% [params_fit2_n(:) sig_fit2_n(1:size(Mxy, 2))] = fit_function(Mxy(:, :,3), TR, flips( :, :), params_fixed_2, params_est_2, [], plot_fits);


			if disp_1xlsx_flag==1
				disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P, R1L, kPL))
				disp(['kPL  = ']); 
				disp('Tumor results:')
				disp([num2str(struct2array(params_fit1_t).')])
				disp('Control results:')
				disp([num2str(struct2array(params_fit1_c).')])
				disp('Noise results:')
				disp([num2str(struct2array(params_fit1_n).')])
			end
			if mod(plot_flag,10)==3
				% figure
					subplot(2,3,3) , plot(t, squeeze(Mxy(1,:,1)),'b')
					hold on , plot(t, squeeze(Mxy(1,:,2)),'r')
					plot(t, squeeze(Mxy(1,:,3)),'g')
					legend('tumor','control','noise')
					title('Pyruvate signals')

					subplot(2,3,6) 
%                     plot(t, squeeze(Mxy(2,:,1)),'b')
					hold on, plot(t, squeeze(Mxy(2,:,2)))
% 					plot(t, squeeze(Mxy(2,:,3)),'g')
					% plot(t, squeeze(sig_fit1_t(:,:)),'bo')
					plot(t, squeeze(sig_fit1_c(:,:)),'o')
					% plot(t, squeeze(sig_fit1_n(:,:)),'go')

					% plot(t, squeeze(sig_fit2_t(:,:)),'b*')
					plot(t, squeeze(sig_fit2_c(:,:)),'*')
					% plot(t, squeeze(sig_fit2_n(:,:)),'g*')
					title(['Lactate signals',filename(1:10)])
                    legend('ori cntl','nml cntl','wo/ns cntl')
			end
		end
			% saveas(fig,['./Andy testing/presentation/Jun 19/',filename(1:6),num2str(A_scaler),' .png']);
	end % loop for one xlsx file
	% save([fitted_dir,tumor_type],'tumor_type','tumor_kpl','errors')
end % loop for one tumor folder






