% Demonstration of single-band hyperpolarized C-13 pulse designs for metabolite specific MRI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2014 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reset SS package globals
%
clear all
ss_opt([]);
ss_globals;
clc

%% lactate-only pulse
fprintf(1, '\nHere''s a C13 single-band excitation pulse to excite only [1-13C]lactate\n');
fprintf(1, 'while not exciting pyruvate, alanine and pyruvate-hydrate, for a 3T clinical system\n\n');

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
          'Spect Correct', 1});
      
% SPECTRAL PULSE PARAMETERS  - large pass/stop bands chosen for wide
% supression regions
B0 = 3e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)
mets(1).name = 'pyr'; 	mets(1).f = -230; 	mets(1).df = 2*df; 		mets(1).ang = 0; 
mets(2).name = 'ala'; 	mets(2).f = -45; 	mets(2).df = 3*df;      mets(2).ang = 0; 
mets(3).name = 'pyrh'; 	mets(3).f = 40; 	mets(3).df = 2*df; 		mets(3).ang = 0; 
mets(4).name = 'lac'; 	mets(4).f = 165; 	mets(4).df = 2*df;      mets(4).ang = 90; 

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 4; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]lac only for 3T clinical system');

fprintf(1,'Hit any key to continue:\n');
pause;

%% lactate-only for high-field, small animal system
fprintf(1, '\nHere is a similar pulse design, but for a 14T small-animal system\n\n');

clear all; ss_opt([]);	ss_globals;			% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 10e-3, ...
	      'Max Grad', 50, ...  % G/cm
	      'Max Slew', 200, ... % G/cm/ms
          'Verse Fraction', 0.5, ...
          'Spect Correct', 1});
      
% SPECTRAL PULSE PARAMETERS  - large pass/stop bands chosen for wide
% supression regions
B0 = 14.1e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)
mets(1).name = 'pyr'; 	mets(1).f = -1080; 	mets(1).df = 2*df; 		mets(1).ang = 0; 
mets(2).name = 'ala'; 	mets(2).f = -210; 	mets(2).df = 2*df;      mets(2).ang = 0; 
mets(3).name = 'pyrh'; 	mets(3).f = 190; 	mets(3).df = 2*df; 		mets(3).ang = 0; 
mets(4).name = 'lac'; 	mets(4).f = 775; 	mets(4).df = 2*df;      mets(4).ang = 90; 

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3.5; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]lac only for 14T animal system');

fprintf(1,'Hit any key to continue:\n');
pause;


%% lactate-only pulse
fprintf(1, '\nFor studies including copolarized 13C-urea or a 13C-urea phantom\n');
fprintf(1, 'Here''s a C13 single-band excitation pulse to excite only [1-13C]lactate\n');
fprintf(1, 'while not exciting pyruvate, alanine and pyruvate-hydrate AND urea, for a 3T clinical system\n\n');

clear all; ss_opt([]);	ss_globals;			% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
          'Spect Correct', 1});
      
% SPECTRAL PULSE PARAMETERS  - large pass/stop bands chosen for wide
% supression regions
B0 = 3e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)
mets(1).name = 'urea'; 	mets(1).f = -470; 	mets(1).df = 2*df;      mets(1).ang = 0; 
mets(2).name = 'pyr'; 	mets(2).f = -230; 	mets(2).df = 2*df; 		mets(2).ang = 0; 
mets(3).name = 'ala'; 	mets(3).f = -45; 	mets(3).df = 3*df;      mets(3).ang = 0; 
mets(4).name = 'pyrh'; 	mets(4).f = 40; 	mets(4).df = 2*df; 		mets(4).ang = 0; 
mets(5).name = 'lac'; 	mets(5).f = 165; 	mets(5).df = 2*df;      mets(5).ang = 90; 

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % minimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = 1;  % thickness (cm)
z_tb = 4; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]lac only for 3T clinical system, avoiding urea');
