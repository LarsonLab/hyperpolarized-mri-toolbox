
clear all;
close all;
clc;

% %%
data_dir = './Andy testing/data/Andy_data/with_vexel_xlsx/';
save_dir = './Andy testing/data/Andy_data/with_vexel_mat/';
[filenames,files_num] = A_get_filenames(data_dir);
disp('filenames: ');
disp(num2str(cell2mat(filenames)));
filenames=char(filenames);

N = 17; % should double check with Renuka Sriram
fit_function = @fit_kPL_perfused_voxel;
line = 4; lac_sheet = 1; flip_sheet = 4;
for file_index=1:files_num
	filename = filenames(file_index,:);
	filename = filename(1:strfind(filename, ' copy.xlsx')-1)
	[flips,Mxy,t,TR,Tin,std_noise,mean_noise,VIF] = A_load_avg_data_perfused([data_dir,filename,' copy.xlsx'],flip_sheet,N);

	save([save_dir,filename],'filename','Mxy','flips','VIF'...
		,'t','TR','Tin','std_noise','mean_noise')

	figure
	subplot(121)
	plot(squeeze(Mxy(1,:,:)))
    title('pyruvate data')
    legend('tumor','control','noise')
    
	subplot(122)
	plot(squeeze(Mxy(2,:,:)))
    title('lactate data')
    legend('tumor','control','noise')
	drawnow, pause(0.5)
end