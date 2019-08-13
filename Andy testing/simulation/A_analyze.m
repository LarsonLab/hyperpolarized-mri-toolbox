clear all; close all; clc;
fitted_dir='./Andy testing/data/perfused_fitted/';
load([fitted_dir,'old_model.mat'])
% old_kpl=tumor_kpl;
clearvars -except fitted_dir old_kpl
load([fitted_dir,'perfused_model_tunning2.mat'])
pf_kpl=tumor_kpl;
[row,col]=size(pf_kpl);
%%
pf_kpl
old_kpl
%%
close all;
a=[1,3,0,0.3];
titles={'no-noi tumor','no-noi cont','no-noi noise', 'noisy tumor','noisy cont','noisy noise'};
fig=figure(5);
set( fig, 'units','normalized', 'outerposition', [0.2 0.1 0.6 0.8], 'Name', 'KPL comparison with outlier');
for i = 1:6
	subplot(2,3,i)
	plot(pf_kpl(:,i),'-*')
	hold on;
	plot(old_kpl(:,i),'-s')
	legend('perfused model','old model')
	% axis(a)
	title(titles{i})
end
% titles={'no-noi tumor','no-noi cont','noi noise', 'noisy tumor','noisy cont','noisy noise'};
%%
fig=figure(6);
a=[1,3,0,1];
set( fig, 'units','normalized', 'outerposition', [0.2 0.1 0.6 0.8], 'Name', 'KPL comparison without outlier');
for i = 1:6
	subplot(2,3,i)
	plot(pf_kpl(1:3,i),'-*')
	hold on;
	plot(old_kpl(1:3,i),'-s')
	legend('perfused model','old model')
    xlabel('xlsx file')
    ylabel('k_{pl}')
	title(titles{i})
	axis(a)
end