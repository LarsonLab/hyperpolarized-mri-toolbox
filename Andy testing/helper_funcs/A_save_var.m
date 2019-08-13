% clear all;
clc;
%%
data_dir = './Andy testing/Andy_data';
tumors={'/786O','/A498','/uok262'};
tumors_num=length(tumors);
tumors=char(tumors);
for tumor_counter=1:tumors_num
    tumor_type=tumors(tumor_counter,:);
end
N = 20;
[filenames,files_num] = A_get_filenames([data_dir,tumor_type]);
filenames=char(filenames);
plot_flag=0;
plot_fits = 0;
plot_mag = 0;
load_avg=1;
%%
size(Mxy)
size(flips)
data_dir = './Andy testing/Andy_data';
tumors={'/786O','/A498','/uok262'};
fitted_dir='./Andy testing/Andy_fitted';

save([fitted_dir,'/',filename],'filename','Mxy','flips','kpl_tmr','kpl_ctr','kpl_n','sig_fit1')
%%

clear all
