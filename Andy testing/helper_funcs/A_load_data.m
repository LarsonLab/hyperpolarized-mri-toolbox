% This file is for Andy to load data from xlsx
function [flips,Mxy,t,TR] = A_load_data(dir,filename,sheet,N)
	% dir = '/Andy testing';
	% filename_test1 = '/book1.xlsx';
	% filename_test2 = '/testing1.xlsx';
	flip_pyr=xlsread([dir,filename],sheet,['S3:S',num2str(3+N-1)])'
	lac_data=xlsread([dir,filename],sheet,['B3:B',num2str(3+N-1)])'
	pyr_data=xlsread([dir,filename],sheet,['C3:C',num2str(3+N-1)])'
	Tin=xlsread([dir,filename],5,'A3')
	TR=xlsread([dir,filename],5,'A4')
	TR=TR-Tin
	% first row pyr
	flips(1, 1: N) = flip_pyr * pi / 180
	% second row lac
	flips(2, 1: N) = ones(1,N) * 90 * pi / 180

	Mxy(1:2, 1:N)=[lac_data;pyr_data]

	t = [0:N - 1] * TR + Tin;