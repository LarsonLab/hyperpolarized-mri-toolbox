% Script to generate a singleband SPSP RF pulse, suitable for imaging
% pyruvate, lactate, and bicarbonate in an HP 1-13C pyruvate experiment


% Set SS package globals
close all, clearvars; clc
ss_globals;
ss_opt([]);				% Reset all options
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 20.5e-3, ...
	      'Num Lobe Iters', 8, ...
	      'Max B1', 1.6, ...
	      'Num Fs Test', 100, ...
	      'Verse Fraction', 0.7, ...
	      'SLR', 0, ...
	      'B1 Verse', 0, ...
	      'Min Order', 0,...
	      'Spect Correct', 1,...
          'Max Grad',5,...
          'Max Slew',20});

%% SPECTRAL PULSE PARAMETERS 
% Design a pulse suitable for imaging C1 pyruvate, lactate, and bicarbonate
% at 3 T
%
% All frequencies are relative to C1 pyruvate
% Note that freq bandwidth is +/-, NOT TOTAL. Need to double to get total bandwidth.
% To generate a pulse suitable for imaging all metabolites, for best
% results place a broad symmetric stopband around the central resonance

disp(['Here''s a C13 single-band excitation pulse for a [1-13C]pyruvate experiment ' ...
    'suitable for imaging C1 pyruvate, lactate, and bicarbonate.  ' ...
    'For best results it places a broad symmetric stopband around the central resonance']);

B0 = 3e4; % G
c13ppm = 1e-6 * B0 * SS_GAMMA; % 1 ppm = gamma_C13 * B0 * 1e-6

% 3T frequency shifts for C1-pyruvate studies
freq_c1 = [395   270   180     0  -322]; %Lac, Pyr-Hydrate, Ala, Pyr, Bic
metabolite_c1 = {'Lactate','Pyruvate-hydrate','Alanine','Pyruvate','Bicarbonate'};
% minimum metabolite separation - using lac to pyr-hyd
fmin_stopband = freq_c1(1) - freq_c1(2);
% maximum metabolite separation, and adding 1 ppm buffer
fmax_stopband = freq_c1(1) - freq_c1(5) + 1*c13ppm;

fcenter_stopband = (fmin_stopband + fmax_stopband)/2;
df_stopband = (fmax_stopband - fmin_stopband)/2;

% metabolite                    frequency (Hz)                  freq bandwidth (Hz)         flip angle (deg)	allowed ripple
mets(1).name = 'passband';      mets(1).f = 0;                  mets(1).df = 1*c13ppm;      mets(1).ang = 90; 	mets(1).d = .02;
mets(2).name = 'stopband1'; 	mets(2).f = fcenter_stopband; 	mets(2).df = df_stopband; 	mets(2).ang = 0; 	mets(2).d = .005;  % low ripple here important to not excite the large pyruvate when exciting other metabolites
mets(3).name = 'stopband2'; 	mets(3).f = -fcenter_stopband;  mets(3).df = df_stopband; 	mets(3).ang = 0; 	mets(3).d = .005; % low ripple here important to not excite the large pyruvate when exciting other metabolites


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
ss_type = 'Flyback Whole'; %spsp type
ptype = 'ex';              %excitation
z_ftype = 'ls';            %spatial filter
z_d1 = 0.01;               %passband ripple
z_d2 = 0.01;               %stopband ripple

% 1) Fs:  868.1 B1: 0.208G Power: 2.831e-05 G^2 ms Dur: 20.6ms
[g,rf,fs] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, [], ss_type, fctr, 0);
%% Simulate SPSP response
% With the SPSP pulse above, confirm frequency selectivity when moving 
% between metabolites. Need to ensure that off-resonance stopband
% excitation <1% since the RF pulse is doing the spectral encoding

close all;
freq_c1 = [395   270   180     0  -322]; %Lac, Pyr-Hydrate, Ala, Pyr, Bic
% Pyr
fplot = freq_c1 - freq_c1(4);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Pyruvate Excitation','NumberTitle','Off')

% Lac
fplot = freq_c1 - freq_c1(1);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Lactate Excitation','NumberTitle','Off')

% Bic
fplot = freq_c1 - freq_c1(5);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','Bicarbonate Excitation','NumberTitle','Off')

% Ala
fplot = freq_c1 - freq_c1(3);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Alanine Excitation','NumberTitle','Off')

fprintf(1,'Hit any key to continue:\n');
pause;

%% Adapt singleband SPSP for alanine
disp(['As seen above, the passband FWHM was too broad to selectively excite alanine. ' ...
    'To create one that is sufficiently narrow, we will need to reduce the passband FWHM ' ...
    'and increase the Max Duration in the SS_OPT structure, ' ...
    'since the total pulsewidth will ultimately limit the passband response ' ...
    '(for Gaussiun modulation, minimum passband FWHM = 2/pw).'])

opt = ss_opt({'Max Duration', 24e-3});
mets(1).df = 0.5*c13ppm;

% Increase stopband width to help suppress ripples from the passband transition

% minimum metabolite separation - pyr-hyd to ala
fmin_stopband = freq_c1(2) - freq_c1(3);
% maximum metabolite separation - lac to bic, and adding 1 ppm buffer
fmax_stopband = freq_c1(1) - freq_c1(5) + 1*c13ppm;

fcenter_stopband = (fmin_stopband + fmax_stopband)/2;
df_stopband = (fmax_stopband - fmin_stopband)/2;

% metabolite                    frequency (Hz)                  freq bandwidth (Hz)         flip angle (deg)	allowed ripple
mets(2).name = 'stopband1'; 	mets(2).f = fcenter_stopband; 	mets(2).df = df_stopband; 	mets(2).ang = 0; 	mets(2).d = .005;  % low ripple here important to not excite the large pyruvate when exciting other metabolites
mets(3).name = 'stopband2'; 	mets(3).f = -fcenter_stopband;  mets(3).df = df_stopband; 	mets(3).ang = 0; 	mets(3).d = .005; % low ripple here important to not excite the large pyruvate when exciting other metabolites

% Regenerate passband & stopband specifications
[fspec, a_angs, d] = create_freq_specs(mets, 0); 

close all;
% 1) Fs:  868.1 B1: 0.134G Power: 1.634e-05 G^2 ms Dur: 24.0ms
[g,rf,fs] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, [], ss_type, fctr, 0);

%% Simulate SPSP response   
% With this updated SPSP pulse, the passband will now be selective enough
% to image alanine as well
close all;

% Pyr
fplot = freq_c1 - freq_c1(4);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Pyruvate Excitation','NumberTitle','Off')

% Lac
fplot = freq_c1 - freq_c1(1);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Lactate Excitation','NumberTitle','Off')

% Bic
fplot = freq_c1 - freq_c1(5);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','Bicarbonate Excitation','NumberTitle','Off')

% Ala
fplot = freq_c1 - freq_c1(3);
ss_plot(g, rf, SS_TS, ptype, z_thk*3, 2.5*[min(fspec) max(fspec)], SS_GAMMA, fplot);
set(gcf,'Name','C1 Alanine Excitation','NumberTitle','Off')


