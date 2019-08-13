% This file is for Andy to load data from xlsx
function [flips,Mxy,t,TR,std_noise] = A_load_vexel_data(dir,line,lac_sheet,flip_sheet,N)
	% dir = '/Andy testing';
	% filename_test1 = '/book1.xlsx';
	% filename_test2 = '/testing1.xlsx';

	lac_data=xlsread(dir,lac_sheet,['D',num2str(line),':',char(double('D')+N-1),num2str(line)]);
	pyr_data=xlsread(dir,lac_sheet+1,['D',num2str(line),':',char(double('D')+N-1),num2str(line)]);
	flip_pyr=xlsread(dir,flip_sheet,['S3:S',num2str(3+N-1)])';
	std_noise=xlsread(dir,sheet,['J',num2str(6+N-1)])';
	
	Tin=xlsread(dir,flip_sheet,'A3');
	TR=xlsread(dir,flip_sheet,'A4');
	TR=TR-Tin;
	% first row pyr
	flips(1, 1: N) = flip_pyr * pi / 180;

	% second row lac
	flips(2, 1: N) = ones(1,N) * 90 * pi / 180;
	flips(2, 1:3) = 0;
	
	Mxy(1:2, 1:N)=[pyr_data;lac_data];

	t = [0:N - 1] * TR + Tin;