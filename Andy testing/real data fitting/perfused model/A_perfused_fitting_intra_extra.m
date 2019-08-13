% Andy's comment:
% 	This file is a modified version of Prof. Larson's test_fit_kPL_fcn.m
%   This file is for testing over the data collected by Renuka Sriram.
% 	Thanks for all the help from Shuyu Tang
% 	Script for testing fit_kPL kinetic model fitting functions

clear all;
close all;
clc;

%%
data_dir = './Andy testing/data/Andy_data/with_vexel_mat/';
fitted_dir='./Andy testing/data/perfused_fitted/';
[filenames,files_num] = A_get_filenames2(data_dir,'.mat');
plot_flag=1;

tumor_kpl=zeros(files_num,6);
para_tumor=zeros(files_num,8); % tumor, control
para_cont=zeros(files_num,8); % tumor, control
% choose fitting function to test
fit_function = @fit_kPL_perfused_voxel;
plot_fits = 0;
%



for file_index=1:files_num % loop for one mat file
	clearvars -except data_dir fitted_dir filenames files_num file_index plot_flag fit_function plot_fits tumor_kpl para_tumor para_cont
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
			legend('tumor','cont','noise','VIF pyr')
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
			legend('tumor','cont','noise','VIF lac')
			% legend('Pyruvate tumor','Lactate tumor','Pyruvate control','Lactate control','Pyruvate noise','Lactate noise')
	end
	std_noise=std_noise*sqrt(pi/2);

	% initial parameter guesses
	R1P_est = 1/25; R1L_est = 1/25; kPL_est = 0.2;  kve_est =0.02; vb_est = 0.2; vif_est=1;
	% kve_est = kve; vb_est = vb;kve = 0.05; vb = 0.1;
	% sum(Mxy(1,:,1))/sum(VIF(1,:))
	disp(['estimate vif from data: ',num2str(sum(Mxy(1,:,1))/sum(VIF(1,:)))]);
	vif_est=sum(Mxy(1,:,1))/sum(VIF(1,:));
	% Test fitting - fit kPL only
	% disp('Fitting kPL, with fixed relaxation rates:')
	% disp('Fixing relaxation rates improves the precision of fitting, but potential')
	% disp('for bias in fits when incorrect relaxation rate is used')
	% disp(' ')

	clear params_fixed params_est params_fit params_fitn_complex params_fitn_mag
	params_fixed.R1P = R1P_est; params_fixed.R1L = R1L_est;
	
	params_est.kPL = kPL_est;
	params_est.S0_P=0.5;
	params_est.S0_L=0.3;
	params_est.VIFscale = vif_est;
	params_est.kve = kve_est; params_est.vb = vb_est;
	% add IV function

	for Dtype = 1:3 % tumor control noise
	    % no noise
	    [params_fit(:,Dtype) Sfit(1:2,1:size(Mxy,2),Dtype), Mxy_ev(1:2,1:size(Mxy,2),Dtype), Mxy_iv(1:2,1:size(Mxy,2),Dtype)] = fit_function(Mxy(:,:,Dtype), VIF, TR, flips(:,:), params_fixed, params_est, [], plot_fits);
	    
	    % % add noise
	    % [params_fitn_complex(:,Dtype) Snfit_complex(1:2,1:size(S,2),  Dtype)] = fit_function(Sn(:,:,Dtype), VIF, TR, flips(:,:,Dtype), params_fixed, params_est, [], plot_fits);

	    % magnitude fitting with noisec
	    % [params_fitn_mag(:,Dtype) Snfit_mag(1:2,1:size(Mxy,2),  Dtype)] = fit_function(Mxy(:,:,Dtype), VIF, TR, flips(:,:),params_fixed, params_est, std_noise, plot_fits);
	end
	%
	flip_descripton{1}='tumor';
	flip_descripton{2}='control';
	flip_descripton{3}='noise';
	description_array = [repmat('    ',3,1),  char(flip_descripton(:))];

	tumor_kpl(file_index,1:3)=[getfield(struct2table(params_fit),'kPL')'];
	para_tumor(file_index,:)=struct2array(params_fit(1));
	para_cont(file_index,:)=struct2array(params_fit(2));
	tumor_kpl

	disp('---------------------------------------------------');
	% disp(sprintf('Input R1 = %f (pyr) %f (lac), kPL = %f', R1P_est, R1L_est, kPL_est))
	disp(filename)
	disp('Noiseless fit results:')
	disp(['kPL  = ']); disp([num2str(getfield(struct2table(params_fit),'kPL')), description_array])
	% disp('Noisy magnitude fit results:')
	% disp(['kPL  = ']); disp([num2str(getfield(struct2table(params_fitn_mag),'kPL')), description_array])
	titles={'tumor','cont', 'noise'};
	% figure
	% plot(t, squeeze(Mxy(1,:,:)))
	% hold on;
	% plot()

if plot_flag==1
    fig2=figure(12);
    set( fig2, 'units','normalized', 'outerposition', [0.2 0.1 0.6 0.8], 'Name', 'Data fitting');
    % set( fig, 'units','normalized', 'outerposition', [0.2 0.1 2 2], 'Name', 'Data fitting');
	subplot(files_num,2,1+(file_index-1)*2)
	plot(t, squeeze(Mxy(1,:,1)))
	hold on;
	plot(t, squeeze(Sfit(1,:,1)),':*')
	% plot(t, squeeze(Snfit_mag(1,:,1)),'--s')
	plot(t, squeeze(Mxy_iv(1,:,1)),':v')
	plot(t, squeeze(Mxy_ev(1,:,1)),':^')
	
	xlabel('time')
	ylabel('pyruvate signal')
	legend('data','fitting','intra vascular','extra vascular')
	title(['file: ',num2str(file_index),', Pyruvate  '])
	subplot(files_num,2,2+(file_index-1)*2) 
	hold on;
	plot(t, squeeze(Mxy(2,:,1)))
	plot(t, squeeze(Sfit(2,:,1)),':*')
	% plot(t, squeeze(Snfit_mag(2,:,1)),'--s')
	plot(t, squeeze(Mxy_iv(2,:,1)),':v')
	plot(t, squeeze(Mxy_ev(2,:,1)),':^')
	
	xlabel('time')
	ylabel('lactate signal')
	legend('data','fitting','intra vascular','extra vascular')
	title(['file: ',num2str(file_index),', Lac fitting'])

end
end
names=fieldnames(params_fit)';
names=names(1:6);
%%
para_tumor=para_tumor(:,1:6);
para_cont=para_cont(:,1:6);
tumor_kpl
load([fitted_dir,'old_model.mat'])
old_kpl
disp(names);
estimate=struct2array(params_est);
disp(['estimate:  ',num2str(estimate)])
para_tumor
para_cont
save([fitted_dir,'perfused_model_tunning2'],'tumor_kpl','names','para_tumor','para_cont','estimate')

%%
figure
hold on;
plot(kPL_est*ones(1,4),'-')
plot(old_kpl(1:4,1),'-o')
plot(tumor_kpl(1:4,1),'-s')
axis([1,4,0,1])
legend('estimate','old','perfused')
xlabel('file index')
ylabel('kpl')
title('Fitted kpl result for real data')