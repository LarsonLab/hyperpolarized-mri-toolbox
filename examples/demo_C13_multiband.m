% Demonstration of multiband hyperpolarized C-13 pulse designs for MRSI/CSI, including explanation
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
clc

%% multiband pulse
fprintf(1, '************************************************************\n')
fprintf(1, 'Here''s a C13 multiband excitation pulse example, for [1-13C]pyr+13C-urea\n');
fprintf(1, 'dynamic MR spectroscopic or chemical shift imaging on a 3T clinical system\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue:\n');
pause;

ss_opt([]);				% Reset all options

% GENERAL PULSE PARAMETERS
ss_type = 'EP Whole';  % Echo-planar design
ptype = 'ex';  % excitation pulse
opt = ss_opt({'Nucleus', 'Carbon', ...
	      'Max Duration', 25e-3, ...
          'Spect Correct', 1});
      
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
s_ftype = 'min';  % minimimum-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+13C-urea Multiband');

fprintf(1,'Hit any key to continue:\n');
pause;

% % for saving pulses:
% 
% ss_save(g,rf,max(a_angs),z_thk, [], 'GE', fspec, a_angs, root_fname);

%% multiband pulse - minimum phase
fprintf(1, '************************************************************\n')
fprintf(1, 'This pulse is shortened by using a Minimum-phase spectral filter\n');
fprintf(1, 'This may add some phase offset between metabolites and may slightly distort lineshapes\n');
fprintf(1, 'but under most circumstances won''t reduce SNR and allows for shorter TEs\n\n');
fprintf(1, 'Here is the resulting pulse with a linear-phase (which is longer and may have higher peak power)\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue:\n');
pause;

s_ftype = 'lin';  % linear-phase spectral filter (will allow for complete refocusing of phase between metabolite peaks with spin-echo)

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+13C-urea Multiband - linear-phase');

fprintf(1,'Hit any key to continue:\n');
pause;

s_ftype = 'min';

%%  minimum phase - spectral correction
fprintf(1, '************************************************************\n')
fprintf(1, 'This pulse also is corrected for chemical-shift slice misregistration\n');
fprintf(1, 'by using the ''Spect Correct'' option.\n');
fprintf(1, 'If the frequency specification bandwidth is too large and/or not sparse, \n');
fprintf(1, 'the spectral correction can fail and/or result in high RF pulse powers.\n\n');
fprintf(1, 'Here is the resulting pulse without spectral correction\n');
fprintf(1, '(the slices for different frequencies are slightly shifted and distorted).\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue:\n');
pause;

opt = ss_opt({'Spect Correct', 0}); % correct for chemical-shift slice misregistration.  This will work for a limited range of frequency shifts 

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+13C-urea Multiband - no spectral correction');

fprintf(1,'Hit any key to continue:\n');
pause;

opt = ss_opt({'Spect Correct', 1}); 
%% multiband pulse - thicker slice
fprintf(1, '************************************************************\n')
fprintf(1, 'Increasing slice thickness from 5mm to 1cm will allow for a higher\n');
fprintf(1, 'spatial time-bandwidth and a sharper slice profile\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue:\n');
pause;

z_thk = 1;  % thickness (cm)
z_tb = 5; % time-bandwidth, proportional to profile sharpness

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+13C-urea Multiband - slab');

fprintf(1,'Hit any key to continue:\n');
pause;

%% short duration pyruvate-lactate pulse

fprintf(1, '************************************************************\n')
fprintf(1, 'Here''s a short C13 multiband excitation pulse, for [1-13C]pyr/lac imaging on a 3T clinical system\n');
fprintf(1, 'with unspecified flip angles for other metabolites for cancer imaging applications.\n');
fprintf(1, 'This type of design is used in UCSF clinical MRSI studies.\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue:\n');
pause;

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
s_ftype = 'min';  % linear-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_angs, d, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);
set(gcf,'Name', '[1-13C]pyr+lac Simple Multiband');


%% short duration pyruvate-lactate pulse with variable flip angles

fprintf(1, '************************************************************\n')
fprintf(1, 'To improve the SNR for dynamic imaging, the flip angle should be varied over time.\n');
fprintf(1, 'Here is an example to design a set of RF pulses, with identical gradients.\n\n');
fprintf(1, 'and different variable flip angles for pyruvate and lactate.\n');
fprintf(1, '************************************************************\n')
fprintf(1,'Hit any key to continue: (design is slower since a set of pulses)\n');
pause;

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
% metabolite			frequency (Hz)		freq bandwidth (Hz)
mets(1).name = 'pyr'; 	mets(1).f = -230; 	mets(1).df = 2*df;
mets(2).name = 'lac'; 	mets(2).f = 165; 	mets(2).df = 2*df; 

% define variable flip angles
% (see hyperpolarized-mri-toolbox for more info and code on schemes)
TR = 3; N = 16; 
kPLest = .05;  R1est = [1/35 1/30]; % [T1pyr, T1lac], 1/s

flips = zeros(2,N);
% pyruvate flips for constant signal
% vfa_const_amp(N, pi/2, exp(-TR * (R1(1)+kPL)))
E1 = exp(-TR * (R1est(1)+kPLest));
flips(1,N) = pi/2;
for n = N-1:-1:1
    flips(1,n) = atan(E1*sin(flips(1,n+1)));
end

% lactate flips for maximum total SNR
% vfa_opt_signal(N, exp(-TR * R1(2)));
E1 = exp(-TR * R1est(2));
flips(2,:) = acos( sqrt((E1^2-E1.^(2.*(N-[1:N]+1))) ./ (1-E1.^(2*(N-[1:N]+1)))) );


dr = .02; % here ripple is a fraction of the flip angle

for n = 1:length(mets)
    fspec(2*n-1) = mets(n).f - mets(n).df;
    fspec(2*n) = mets(n).f + mets(n).df;
end
fmid = (mets(1).f+mets(end).f)/2;
fspec = fspec - fmid;

fctr = 0;  % force pulse design to optimize for center of frequency specification
s_ftype = 'min';  % linear-phase spectral filter

% SPATIAL PULSE PARAMETERS
z_thk = .5;  % thickness (cm)
z_tb = 3; % time-bandwidth, proportional to profile sharpness
z_ftype='ls';  % least-squares filter design
z_d1 = 0.01;  z_d2 = 0.01;  % slice profile pass and stop-band ripples, respectively

% DESIGN THE PULSE!
[g,rf,fs,z,f,mxy] = ...
    ss_design_dyn(z_thk, z_tb, [z_d1 z_d2], fspec, flips, dr, ptype, ...
	      z_ftype, s_ftype, ss_type, fctr);

% % for saving dynamic pulses:
% for t = 1:N
%     ss_save_dyn(g(:,t),rf(:,t),max(flips(:,t)),z_thk, [], 'GE', fspec, flips(:,t), root_fname, t);
% end


