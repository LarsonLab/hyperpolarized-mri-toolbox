% Script to test slice and frequency shifting
% Uses a singleband SPSP RF pulse, suitable for imaging
% pyruvate, lactate, and bicarbonate in an HP 1-13C pyruvate experiment

close all, clearvars; clc
ss_globals;
ss_opt([]);				% Reset all options
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 30e-3, ...
	      'Num Lobe Iters', 10, ...
	      'Max B1', 1.6, ...
	      'Num Fs Test', 100, ...
	      'Verse Fraction', 0.7, ...
	      'SLR', 0, ...
	      'B1 Verse', 0, ...
	      'Min Order', 1,...
	      'Spect Correct', 1,...
          'Max Grad',5,...
          'Max Slew',20});

% SPECTRAL PULSE PARAMETERS 
% Design a pulse suitable for imaging C1 pyruvate, lactate, and bicarbonate
% at 3 T
%
% All frequencies are relative to C1 pyruvate
% Note that freq bandwidth is +/-, NOT TOTAL. Need to double to get total bandwidth.
% To generate a pulse suitable for imaging all metabolites, for best
% results place a broad symmetric stopband around the central resonance
B0 = 3e4; % G
c13ppm = 1e-6 * B0 * SS_GAMMA; % 1 ppm = gamma_C13 * B0 * 1e-6

% 3T frequency shifts for C1-pyruvate studies
freq_c1 = [395   270   180     0  -322]; %Lac, Pyr-Hydrate, Ala, Pyr, Bic
metabolite_c1 = {'Lactate','Pyruvate-hydrate','Alanine','Pyruvate','Bicarbonate'};

% metabolite                    frequency (Hz)          freq bandwidth (Hz)         flip angle (deg)	allowed ripple
mets(1).name = 'passband';      mets(1).f = 0;          mets(1).df = 1*c13ppm;      mets(1).ang = 90; 	mets(1).d = .01;

fmin_stopband = freq_c1(1) - freq_c1(2);
fmax_stopband = freq_c1(1) - freq_c1(5) + 1*c13ppm;
fcenter_stopband = (fmin_stopband + fmax_stopband)/2;
df_stopband = (fmax_stopband - fmin_stopband)/2;

mets(2).name = 'lac_bic_stopband1'; 	mets(2).f = fcenter_stopband; 	mets(2).df = df_stopband; 		mets(2).ang = 0; 	mets(2).d = .005;  % low ripple here important to not excite the large pyruvate when exciting other metabolites
mets(3).name = 'bic_lac_stopband2'; 	mets(3).f = -fcenter_stopband;  mets(3).df = df_stopband; 		mets(3).ang = 0; 	mets(3).d = .005; % low ripple here important to not excite the large pyruvate when exciting other metabolites

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets, 0);
fctr = 0;         % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
% The SPSP type has the biggest impact on the overall response.
% Echoplanar designs will have a larger Fs and can provide thinner slices,
% but they are more sensitive to RF/gradient timing errors and the opposed
% null at Fs/2 can make it difficult to find a pulse suitable for imaging
% all metabolites
z_thk = 1.25;              %slice thickness (cm)
z_tb = 2;                  %spatial time bandwidth
ss_type = 'Flyback Whole';      %spsp type
% ss_type = 'EP Whole';      %spsp type - this should also work for echo-planar design
ptype = 'ex';              %excitation
z_ftype = 'ls';            %spatial filter
z_d1 = 0.01;               %passband ripple
z_d2 = 0.01;               %stopband ripple

% 1) Fs:  868.1 B1: 0.223G Power: 3.051e-05 G^2 ms Dur: 20.6m
[g,rf,fs] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, [], ss_type, fctr, 0);

%% test frequency shifing

close all;

for Ifreq = 1:length(freq_c1)
    rf_shift = ss_shift(g,rf,0,freq_c1(Ifreq));
    ss_plot(g, rf_shift, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, freq_c1);
    set(gcf,'Name',[metabolite_c1{Ifreq} ' Excitation'],'NumberTitle','Off')

end

%% test slice shifting
num_slices = 5;
slice_gap = 0;
z_shift = (z_thk+slice_gap)*([1:num_slices] - (num_slices+1)/2);

fplot = freq_c1 - freq_c1(4);
for Islice = 1:num_slices
    rf_shift = ss_shift(g, rf, z_shift(Islice));
    ss_plot(g, rf_shift, SS_TS, ptype, z_thk*num_slices, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
end

