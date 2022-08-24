%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2011 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reset SS package globals
%
ss_opt([]);
ss_globals;


%%
clc
fprintf(1, '\nHere''s a multiband spin-echo pulse for MR spectroscopy at 1.5T\n');
fprintf(1, 'It includes fat suppression and partial water suppression\n');
fprintf(1, 'The specification includes: \n');
fprintf(1, '\t-no refocusing of fat (0-degree flip at ~ -220 Hz)\n');
fprintf(1, '\t-full refocusing of choline and citrate (180-degree flip from -75 to -140 Hz)\n');
fprintf(1, '\t-partial refocusing of water (~37-degree flip at 0 Hz)\n');
fprintf(1, '(similar to Schricker et al, Magn Reson Med 46: 1079-1087 (2001), DOI: 10.1002/mrm.1302 )\n');
fprintf(1,'Hit any key to continue:\n');
pause

% Water/fat chemical shifts (ppm)
%
df = 0.5e-6;				% Conservative +-0.5 ppm shim requirement

fat1 = 0.9e-6;
fat2 = 1.3e-6;
cho = 3e-6;
water = 4.7e-6;

% Convert to frequency specification
%
B0 = 15000;
gamma_h = 4258;
% bands: fat, cho, water
fspec = B0 * ([(fat1-df) (fat2+df) (cho-df) (cho+df) (water-df) (water+df)]-water) * gamma_h;
mxy = [0 1 0.1];  % refocused mxy for a single 180 pulse
d = [0.005 0.01 0.005]; % mxy ripple values

fat_ctr = ((fat1+fat2)/2-water)*B0 * gamma_h;
cho_ctr = (cho-water)*B0 * gamma_h;
water_ctr = 0;

% Set up pulse parameters
%
z_thk = 1; % cm
z_tb = 5; % time-bandwidth

% Set up spectral/spatial specifications
%

a_ang = 2*asin(sqrt(mxy));  % convert Mxy amplitude to flip angle
ptype = 'se';  % pulse type: 'se' spin-echo
z_ftype='ls';				% Use this to get rid of "Conolly Wings" 
z_d1 = 0.01;
z_d2 = 0.01;
f_ctr = cho_ctr;

s_ftype = 'max';			% max-phase spectral 
ss_type = 'Flyback Whole';
dbg = 0;				% dbg level: 0 -none, 1 - little, 2 -lots, ...

ss_opt([]);
opt = ss_opt({'Max Duration', 30e-3, ...  % ms
	      'Num Lobe Iters', 5, ...
	      'Spect Correct', 0, ...
	      'SLR', 1, ...
	      'Verse Fraction', 0.8, ...
	      'Max Grad', 4, ...  % these gradient specifications could be updated for modern systems
	      'Max Slew', 15, ...
          'Min Order', 1})


[g,rf,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'MRS Dual-band spin-echo pulse');

fprintf(1,'Hit any key to continue:\n');
pause

%%
fprintf(1, '\nBy default the spectral-spatial package uses the shortest possible RF pulse,\n');
fprintf(1, 'but this can cause large ripples outside of the specified bands, as seen in this pulse.\n');
fprintf(1, 'The ripples, and also the peak and total power, are reduced by turning off the\n');
fprintf(1, 'Min Order option, in which case the resulting pulses are the Max Duration:\n');

opt = ss_opt({'Min Order', 0})
[g,rf,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a_ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'MRS Dual-band spin-echo pulse - Reduced power, longer duration');
