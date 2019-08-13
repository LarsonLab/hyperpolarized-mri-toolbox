% Andy's comment:
% 	This file is a modified version of Prof. Larson's test_fit_kPL_fcn.m
%   This file is for testing over the data collected by Renuka Sriram.
% 	Thanks for all the help from Shuyu Tang
% 	Script for testing fit_kPL kinetic model fitting functions

clear all;
close all;
clc;

%%
data_dir = './Andy testing/others/Andy_data/with_vexel_mat/';
fitted_dir='./Andy testing/others/perfused_fitted/';
[filenames,files_num] = A_get_filenames2(data_dir,'.mat');
plot_flag=1;

% tumor_kpl=zeros(files_num,2,3);
fit_parameter=zeros(files_num,2,3);
% fit_signal=zeros(files_num,2,3);
% para_cont=zeros(files_num,2,3); % tumor, control
% choose fitting function to test
fit_function = @changed_fit_kPL_perfusion;
plot_fits = 0;
%



for file_index=1:files_num % loop for one mat file
	clearvars -except data_dir fitted_dir filenames files_num file_index plot_flag fit_function plot_fits tumor_kpl fit_parameter para_cont
	% file_index=1; % for debug
	load([char(filenames(file_index))])
	if plot_flag==1
		fig1=figure(11);
		    set( fig1, 'units','normalized', 'outerposition', [0.2 0.1 0.6 0.8], 'Name', 'Data fitting');

			% plot pyr and lac
			subplot(files_num,2,file_index*2-1) , plot(t, squeeze(Mxy(1, : , 1)))
			hold on
			plot(t, squeeze(Mxy(1, :,2 )),':')
			plot(t, squeeze(Mxy(1, :,3 )),'-')
			plot(t, squeeze(VIF(1,:)),'-*')
			legend('tumor','cont','noise','VIF pyr','Location','northeastoutside')
			% hold off;

	xlabel('time')
	ylabel('pyruvate signal')
			title(['file: ',num2str(file_index),', Pyruvate  ',filename(1:10),' data'])
			

			subplot(files_num,2,file_index*2); hold on;
			plot(t, squeeze(Mxy(2, :,1 )))
			plot(t, squeeze(Mxy(2, :,2 )),':')
			plot(t, squeeze(Mxy(2, :,3 )),'-')
			plot(t, squeeze(VIF(2,:)),'-*')
			% hold off;

	xlabel('time')
	ylabel('lactate signal')
			title(['file: ',num2str(file_index),', Lactate data'])
			legend('tumor','cont','noise','VIF lac','Location','northeastoutside')
			% legend('Pyruvate tumor','Lactate tumor','Pyruvate control','Lactate control','Pyruvate noise','Lactate noise')
	end
	std_noise=std_noise*sqrt(pi/2);

	% initial parameter guesses
	R1P_est = 1/25; R1L_est = 1/25; kPL_est = 0.2;  kve_est =0.02; vb_est = 0.2; vif_est=1;
	% kve_est = kve; vb_est = vb;kve = 0.05; vb = 0.1;
	% sum(Mxy(1,:,1))/sum(VIF(1,:))
	disp(['estimate vif from data: ',num2str(sum(Mxy(1,:,1))/sum(VIF(1,:)))]);
	vif_est=sum(Mxy(1,:,1))/sum(VIF(1,:))
	% Test fitting - fit kPL only
	% disp('Fitting kPL, with fixed relaxation rates:')
	% disp('Fixing relaxation rates improves the precision of fitting, but potential')
	% disp('for bias in fits when incorrect relaxation rate is used')
	% disp(' ')

	clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
	params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
	
	params_est.kPL = kPL_est;
	% params_est.S0_P=0.5;
	% params_est.S0_L=0.3;
	params_est.VIFscale = vif_est;
	params_fixed.kve = kve_est;
	 % params_est.vb = vb_est;
	% add IV function

	for Dtype = 1:3 % tumor control noise
		% Dtype=1
	    % no noise
	    % [params_fit(:,Dtype) Sfit(1:2,1:size(Mxy,2),Dtype), Mxy_ev(1:2,1:size(Mxy,2),Dtype), Mxy_iv(1:2,1:size(Mxy,2),Dtype)] = fit_function(Mxy(:,:,Dtype), VIF(1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
	    [params_fit(:,Dtype) Sfit(:,Dtype)] = fit_function((Mxy(:,:,Dtype)), VIF(1,:), TR, flips(:,:), params_fixed, params_est, [], plot_fits);
	    
	    % magnitude fitting with noisec
	    % [params_fitn_mag(:,Dtype) Snfit_mag(1:2,1:size(Mxy,2),  Dtype)] = fit_function(Mxy(:,:,Dtype), VIF, TR, flips(:,:),params_fixed, params_est, std_noise, plot_fits);
	end
	fit_parameter(file_index,:,:)=[getfield(struct2table(params_fit),'kPL')';getfield(struct2table(params_fit),'VIFscale')'];
	

if plot_flag==1
    fig2=figure(12);
    set( fig2, 'units','normalized', 'outerposition', [0.2 0.1 0.6 0.8], 'Name', 'Data fitting');
	subplot(files_num,2,1+(file_index-1)*2)
	plot(t, squeeze(Mxy(1,:,:)))
	xlabel('time')
	ylabel('pyruvate signal')
	legend('tumor data','control data','noise data','Location','northeastoutside')
	title(['file: ',num2str(file_index),', Pyruvate  ',filename(1:10),' fitting'])
	subplot(files_num,2,2+(file_index-1)*2) 
	hold on;
	plot(t, squeeze(Mxy(2,:,:)),'--o')
	plot(t, squeeze(Sfit(:,:)),':*')
	xlabel('time')
	ylabel('lactate signal')
	legend('tumor data','control data','noise data','tumor fitting','control fitting','noise fitting','Location','northeastoutside')
	title(['file: ',num2str(file_index),', Lac fitting'])


	figure(13)
	subplot(1,2,1)
	plot(squeeze(fit_parameter(:,1,:)))
	legend('tumor','contorl','noise')
	title('kpl')
	subplot(1,2,2)
	plot(squeeze(fit_parameter(:,2,:)))
	legend('tumor','contorl','noise')
	title('VIFscale')
end
end