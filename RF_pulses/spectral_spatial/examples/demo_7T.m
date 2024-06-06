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


%% Reset SS package globals
%

ss_globals;
fprintf(1, '************************************************************\n')
fprintf(1, 'Here are some parameters that can be set:\n');
fprintf(1, 'You can type "help ss_globals" for more information on\n');
fprintf(1, 'each of them.\n');
fprintf(1, '************************************************************\n')
fprintf(1, '\n');
ss_opt([])


%% Set up small-tip water/fat spectral/spatial
%

% Water/fat chemical shifts
%
df = 0.5e-6;				% Conservative shim requirement
water = 4.7e-6;			
fat2 = 1.3e-6;
fat1 = 0.9e-6;

% Convert to frequency
%
B0 = 70000;
gamma_h = 4258;
fspec = B0 * ([(fat1-df) (fat2+df) (water-df) (water+df)]-water) * gamma_h;

water_ctr = (fspec(3) + fspec(4))/2;
fat_ctr = (fspec(1) + fspec(2))/2;


% Set up pulse parameters
%
ang = pi/6;
z_thk = 1;
z_tb = 4;

% Set up spectral/spatial specifications
%
a = [0 1];
d = [0.01 0.01];
ptype = 'ex';
z_ftype='ls';				% Use this to get rid of "Conolly Wings" 
z_d1 = 0.01;
z_d2 = 0.01;
f_ctr = [];

s_ftype = 'min';			% min-phase spectral 
ss_type = 'EP Whole';
dbg = 0;				% dbg level
                                        % 0 -none, 1 - little, 2 -lots, ...

default_opt = {'Max Duration', 20e-3, ...
          'Max Grad', 4, 'Max Slew', 20, ...  % 40 mT/m, 200 mT/m/ms
	      'Num Lobe Iters', 5, ...
          'Min Order', 1, ...
	      'Spect Correct', 0, ...
	      'SLR', 0, ...
	      'Verse Fraction', 0.9};  

opt = ss_opt(default_opt);
      
% spectral correction
ss_opt({'Spect Correct', 1, ...
    'Spect Correct Reg', 0.001});  % small amount of regularization added to better condition spectral correction inversion



fprintf(1, '\n************************************************************\n')
fprintf(1, 'Here''s an example of a water/fat spectral spatial pulse for 7T\n\n');
fprintf(1, 'The pulse is a echo-planar (''EP'') design to accomodate the large\n');
fprintf(1, 'spectral bandwidth at ultra-high field.\n');
fprintf(1, 'However, EP pulses are more sensitive to eddy currents and timing errors.\n');
fprintf(1, 'The frequency spec includes a passband at [%3f,%3f]\n',fspec(3),fspec(4));
fprintf(1, 'and stopbands at [%3f,%3f]\n',fspec(1),fspec(2));
fprintf(1, 'This design also reduces stopband ripple using the ''Spect Correct'' option.\n');
fprintf(1, 'Here is an example of a fat-water EP trajectory design\n');
fprintf(1, '************************************************************\n');

[g_ew,rf_ew,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);

set(gcf,'Name', 'Water/Fat EP Whole - Spect Correct');




