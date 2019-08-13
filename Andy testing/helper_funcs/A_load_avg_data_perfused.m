% This file is for Andy to load data from xlsx
function [flips,Mxy,t,TR,Tin,std_noise,mean_noise,VIF] = A_load_avg_data_perfused(dir,sheet,N)

	lac_data_tumor=xlsread(dir,sheet,['B6:B',num2str(6+N-1)])';
	pyr_data_tumor=xlsread(dir,sheet,['C6:C',num2str(6+N-1)])';
	lac_data_control=xlsread(dir,sheet,['F6:F',num2str(6+N-1)])';
	pyr_data_control=xlsread(dir,sheet,['G6:G',num2str(6+N-1)])';
	lac_data_noise=xlsread(dir,sheet,['N6:N',num2str(6+N-1)])';
	pyr_data_noise=xlsread(dir,sheet,['O6:O',num2str(6+N-1)])';

	mean_noise=xlsread(dir,sheet,['N',num2str(5+N-1+3)])';
	std_noise=xlsread(dir,sheet,['N',num2str(6+N-1+3)])';

	VIF(1,:)=xlsread(dir,sheet,['J6:J',num2str(6+N-1)])';
	VIF(2,:)=xlsread(dir,sheet,['K6:K',num2str(6+N-1)])';

	flip_pyr=xlsread(dir,sheet,['S6:S',num2str(6+N-1)])';
	Tin=xlsread(dir,sheet,'A3');
	TR=xlsread(dir,sheet,'A4');
	TR=TR-Tin;
	% first row pyr
	flips(1, 1: N) = flip_pyr * pi / 180;

	% second row lac
	flips(2, 1: N) = ones(1,N) * 90 * pi / 180;
	% flips(2, 1:) = 0;
	Mxy(1:2, 1:N,1)=[pyr_data_tumor;lac_data_tumor];
	Mxy(1:2, 1:N,2)=[pyr_data_control;lac_data_control];
	Mxy(1:2, 1:N,3)=[pyr_data_noise;lac_data_noise];
	t = [0:N - 1] * TR + Tin;