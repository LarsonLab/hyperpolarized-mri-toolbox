% Demonstration of hyperpolarized C-13 pulse designs, including explanation
% of some key pulse characteristics

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


%% multiband pulse
fprintf(1, '\nHere''s a C13 multiband excitation pulse, for [1-13C]pyr+13C-urea dynamic imaging on a 3T clinical system\n\n');

ss_opt([]);				% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';  % Echo-planar design
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
	      'Max B1', 0.5});
      
% SPECTRAL PULSE PARAMETERS 
B0 = 3e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)	allowed ripple
mets(1).name = 'urea'; 	mets(1).f = -465; 	mets(1).df = 1.5*df; 		mets(1).ang = 6; 	mets(1).d = .005;
mets(2).name = 'pyr'; 	mets(2).f = -230; 	mets(2).df = 2*df; 		mets(2).ang = 6; 	mets(2).d = .002;
mets(3).name = 'ala'; 	mets(3).f = -45; 	mets(3).df = 1.5*df; 	mets(3).ang = 12; 	mets(3).d = .005;
mets(4).name = 'pyrh'; 	mets(4).f = 40; 	mets(4).df = 1*df; 		mets(4).ang = 2.5; 	mets(4).d = .015;
mets(5).name = 'lac'; 	mets(5).f = 165; 	mets(5).df = 1.5*df; 	mets(5).ang = 12; 	mets(5).d = .005;

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'lin';  % linear-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

fprintf(1,'Hit any key to continue:\n');
pause;

%% multiband pulse - minimum phase
fprintf(1, '\nThe pulse duration can be shortened by using a Minimum-phase spectral filter\n');
fprintf(1, 'This will add some phase offset between metabolites and may slightly distort lineshapes\n');
fprintf(1, 'But under most circumstances won''t reduce SNR and allows for shorter TEs\n\n');

s_ftype = 'min';  % minimum-phase spectral filter (will add phase between metabolite peaks, use 'lin' if this is undesireable)

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

fprintf(1,'Hit any key to continue:\n');
pause;

%% multiband pulse - thicker slice
fprintf(1, '\nIncreasing slice thickness from 5mm to 1cm will allow for a higher\n');
fprintf(1, 'spatial time-bandwidth and a sharper slice profile\n\n');

z_thk = 1;  % thickness (cm)
z_tb = 5; % time-bandwidth, proportional to profile sharpness

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

fprintf(1,'Hit any key to continue:\n');
pause;
%%  minimum phase - spectral correction
fprintf(1, '\nTurning on the ''Spect Correct'' option will to a limited extent correct\n');
fprintf(1, 'chemical-shift slice misregistration.\n');
fprintf(1, 'NOTE: If the frequency specification bandwidth is too large, \n');
fprintf(1, 'this can fail and result in high RF pulse powers.\n\n');

opt = ss_opt({'Spect Correct', 1}); % correct for chemical-shift slice misregistration.  This will work for a limited range of frequency shifts 

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

fprintf(1,'Hit any key to continue:\n');
pause;

%% short duration pyruvate-lactate pulse

fprintf(1, '\nHere''s a short C13 multiband excitation pulse, for [1-13C]pyr/lac imaging on a 3T clinical system\n');
fprintf(1, 'with unspecified flip angles for other metabolites.  Similar to a design used in first-in-man clinical trial.\n\n');

clear all; ss_opt([]); ss_globals;				% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 8e-3, ...
	      'Max B1', 0.5, ...
          'Spect Correct', 1});
      
% SPECTRAL PULSE PARAMETERS 
B0 = 3e4; % G
df = 0.5e-6 * B0 * SS_GAMMA; % 0.5 ppm = gamma_C13 * B0 * 0.5e-6
% metabolite			frequency (Hz)		freq bandwidth (Hz)		flip angle (deg)	allowed ripple
mets(1).name = 'pyr'; 	mets(1).f = -230; 	mets(1).df = 2*df; 		mets(1).ang = 6; 	mets(1).d = .005;
mets(2).name = 'lac'; 	mets(2).f = 165; 	mets(2).df = 2*df;      mets(2).ang = 12; 	mets(2).d = .005;

% create vectors of angles, ripples, and band edges for input to pulse design
[fspec, a_angs, d] = create_freq_specs(mets);
fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'lin';  % linear-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

fprintf(1,'Hit any key to continue:\n');
pause;


%% lactate-only pulse
fprintf(1, '\nHere''s a C13 single-band excitation pulse, for a 3T clinical system\n\n');

clear all; ss_opt([]);	ss_globals;			% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
	      'Max B1', 0.5, ...
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

fprintf(1,'Hit any key to continue:\n');
pause;

%% lactate-only for high-field, small animal system
fprintf(1, '\nHere''s a C13 single-band excitation pulse, for a 14T small-animal system\n\n');

clear all; ss_opt([]);	ss_globals;			% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 10e-3, ...
	      'Max B1', 1.5, ... % G
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
