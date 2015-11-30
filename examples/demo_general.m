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
fprintf(1, 'Hit key to continue:');
pause;


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
B0 = 15000;
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
d = [0.02 0.005];
ptype = 'ex';
z_ftype='ls';				% Use this to get rid of "Conolly Wings" 
z_d1 = 0.01;
z_d2 = 0.01;
f_ctr = [];

s_ftype = 'min';			% min-phase spectral 
ss_type = 'Flyback Half'; 	        % Flyback, symmetric frequency
dbg = 0;				% dbg level
                                        % 0 -none, 1 - little, 2 -lots, ...

default_opt = {'Max Duration', 16e-3, ...
	      'Num Lobe Iters', 5, ...
          'Min Order', 1, ...
	      'Spect Correct', 0, ...
	      'SLR', 0, ...
	      'Verse Fraction', 0.9};

opt = ss_opt(default_opt);
      
    
fprintf(1, '\n************************************************************\n')
fprintf(1, 'Here''s an example of a water/fat spectral spatial pulse for 1.5T\n\n');
fprintf(1, 'The pulse is a ''Flyback Half'' type, meaning that it\n');
fprintf(1, 'uses a flyback trajectory, and has a asymmetric frequency\n');
fprintf(1, 'response.\n');
fprintf(1, 'The frequency spec includes a passband at [%3f,%3f]\n',fspec(3),fspec(4));
fprintf(1, 'and stopbands at [%3f,%3f]\n',fspec(1),fspec(2));
fprintf(1, '************************************************************\n')
fprintf(1, '\n');

[g,rf,fs,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);

set(gcf,'Name', 'Water/Fat Flyback Half');

fprintf(1,'Hit any key to continue:\n');
pause;
%%
fprintf(1, '\n************************************************************\n')
fprintf(1, 'In contrast, here''s a Flyback Half pulse with symmetric response\n');
fprintf(1, 'but with the response specified to be centered on water.\n');
fprintf(1, 'Note that the pulse has a much different response from the\n');
fprintf(1, 'previous design in which the center frequency was not specified.\n');
fprintf(1, '************************************************************\n')
fprintf(1, '\n');

f_ctr = water_ctr;
[g_water,rf_water,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat Flyback Half - Fat Center');

fprintf(1,'Hit any key to continue:\n');
pause;

f_ctr = [];
%%

fprintf(1, '\n************************************************************\n')
fprintf(1, 'We can also emulate a conventional low-pass spectral filter design\n');
fprintf(1, 'as a comparison to our multiband spectral specification\n');
fprintf(1, 'Note that the transition region on the multiband is narrower\n');

f_ctr = water_ctr;

% get sampling freq from previous soln and duration of rf
%
zidx = find(g_water == 0);
fs_max = 1 / (zidx(2) * SS_TS);
flp = fspec;
flp(1) = -0.5*fs_max + f_ctr;
ss_opt({'Max Duration', 11e-3});


[g_lp,rf_lp,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], flp, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat Flyback Half - LP Emulate');

fprintf(1, 'The power is lower too: %f compared to %f\n', ...
	sum(abs(rf).^2), sum(abs(rf_lp).^2));
fprintf(1, '************************************************************\n')

% Plot spectral responses
%
z = 0;
f = linspace(-fs_max/2, fs_max/2, 200);
mxy = ss_response_mxy(g,rf,z,f,SS_TS,SS_GAMMA);
mxy_lp = ss_response_mxy(g_lp,rf_lp,z,f,SS_TS,SS_GAMMA);

figure;
plot([f(:),f(:)],abs([mxy(:) mxy_lp(:)]));
hold on;
a2 = [a(:).'; a(:).'];
a2 = a2(:).';
plot_spec(fspec,sin(ang)*a2,d);
title('Spectral Response');
xlabel('Frequency');
ylabel('M_{xy}');
legend('Multiband', 'LP');

fprintf(1,'Hit any key to continue:\n');
pause;

%%

fprintf(1, '\n************************************************************\n')
fprintf(1, 'By default the spectral-spatial package uses the shortest possible\n');
fprintf(1, 'RF pulse, but this can cause large ripples outside of the specified \n');
fprintf(1, 'bands which also results in higher peak and total power.\n');
fprintf(1, 'The ripples and power can be reduced by turning off the\n');
fprintf(1, '''Min Order'' option, in which case the resulting pulses will become longer\n');
fprintf(1, '(approximately the ''Max Duration'')\n');
fprintf(1, '************************************************************\n')
fprintf(1, '\n');

opt = ss_opt({'Min Order', 0});
[g,rf,fs,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);

set(gcf,'Name', 'Water/Fat Flyback Half - Longer duration');

fprintf(1,'Hit any key to continue:\n');
pause;
%%

fprintf(1, '\n************************************************************\n')
fprintf(1, 'The pulse design will fail if the specifications are too aggressive\n');
fprintf(1, 'In the following case, the Max Duration is too short\n');
fprintf(1, 'ss_design returns an approximate solution with increased ripples.\n');
fprintf(1, '(This could be corrected by reducing bandwidths, increasing ripple,\n');
fprintf(1, 'or increasing the Max Duration.)\n');
fprintf(1, '************************************************************\n')
fprintf(1, '\n');
fprintf(1,'Hit any key to continue:\n');
pause;

opt = ss_opt({'Max Duration', 6e-3, ...
    'Min Order', 1});
[g,rf,fs,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);

set(gcf,'Name', 'Water/Fat Flyback Half - Increased ripple');

fprintf(1,'Hit any key to continue:\n');
pause;
%%

fprintf(1, '\n************************************************************\n')
fprintf(1, 'Here''s an example of how we can verse the gradient lobes\n');
fprintf(1, 'to reduce the peak B1\n');
fprintf(1, '************************************************************\n')
fprintf(1, '\n');

opt = ss_opt(default_opt);
opt = ss_opt({'B1 Verse', 1, ...
	      'Max B1', 0.02});

[g,rf,fs,z,f,mxy] = ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
			    z_ftype, s_ftype, ss_type, f_ctr, dbg);

set(gcf,'Name', 'Water/Fat Flyback Half - B1 Min');

fprintf(1,'Hit any key to continue:\n');
pause;

opt = ss_opt({'B1 Verse', 0});

%%

fprintf(1, '\n************************************************************\n')
fprintf(1, 'What happens if we do not require a symmetric frequency\n');
fprintf(1, 'response?  Good question.  We get a shorter pulse still!!\n');
fprintf(1, '************************************************************\n')

ss_opt(default_opt);
ss_type = 'Flyback Whole';
f_ctr = [];

[g_fw,rf_fw,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat Flyback Whole');

fprintf(1,'Hit any key to continue:\n');
pause;

%%
fprintf(1, '\n************************************************************\n');
fprintf(1, 'The package also allows for echo-planar (''EP'') trajectory designs\n');
fprintf(1, 'which can reduce the pulse duration and allow for thinner slices\n');
fprintf(1, 'but are more sensitive to eddy currents and timing errors.\n');
fprintf(1, 'Here is an example of the original design using a EP trajectory\n');
fprintf(1, 'There is some more ripple in the stopband, but it still meets the spec\n');
fprintf(1, '************************************************************\n');

ss_opt(default_opt);
%opt = ss_opt({'Max B1', 0.4});
ss_type = 'EP Whole';
f_ctr = [];
z_thk = .5;

%d = [0.01 0.005];

[g_ew,rf_ew,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat EP Whole');

fprintf(1,'Hit any key to continue:\n');
pause;

%%

fprintf(1, '\n************************************************************\n');
fprintf(1, 'We can design the same pulse centered on the stopband\n');
fprintf(1, 'to reduce the ripple, but this won''t work well when\n');
fprintf(1, 'we have multiple stopbands\n');
fprintf(1, '************************************************************\n');


ss_type = 'EP Whole';
f_ctr = fat_ctr;
[g_fw,rf_fw,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat EP Whole - Fat Ctr');

fprintf(1,'Hit any key to continue:\n');
pause;



%%

fprintf(1, '\n************************************************************\n');
fprintf(1, 'There is a more sophisticated solution....\n');
fprintf(1, 'The ripple arises from the nonuniform sampling in time\n');
fprintf(1, 'that occurs for EP trajectories\n');
fprintf(1, '************************************************************\n');

ss_type = 'EP Whole';
f_ctr = [];
ss_opt({'Spect Correct', 1});
[g_fw,rf_fw,fs,z,f,mxy] = ...
    ss_design(z_thk, z_tb, [z_d1 z_d2], fspec, a*ang, d, ptype, ...
	      z_ftype, s_ftype, ss_type, f_ctr, dbg);
set(gcf,'Name', 'Water/Fat EP Whole - Spect Correct');

fprintf(1,'Hit any key to continue:\n');
pause;

