%%%% Example: 2 site pyr->lac, pyruvate: 2D GRE, lactate: 3D bSSFP  %%%%%%

%set parameters
clear; close all;

%flip angles
flips.P = 20;
FAL = 60;
Nz = 16; interleaves = 4;
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60];

% TR
TRP = .18; TRL = .01529; %s 
TR.P = TRP * 16;
TR.L = repmat(TRL, [1, Nz*interleaves]);
cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL];

%relaxation
R1 = [1/30, 1/25];
R2.L = 1/1;
k = [0.002];

TempRes =4;
Nt= 30;
Mz0 = [1 0.1];
spoilers = [1 0];
plot_flag = 1;


%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'spoilers', spoilers, 'plot_flag', plot_flag);

%% %%%% Example: 2 site pyr->lac, pyruvate: 3D GRE, lactate: 3D bSSFP  %%%%%%

%set parameters
clear; close all;

%flip angles
FAP = 3;
FAL = 60;
Nz = 16; interleaves = 4;
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60];

% TR
TRP = .18; TRL = .01529; %s 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz*interleaves]);
cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL];

%relaxation
R1 = [1/30, 1/25];
R2.L = 1/1;
k = [0.002]; %[0.02];

TempRes = 4;
Nt= 30;
Mz0 = [1 0.1];
spoilers = [1 0];
plot_flag = 1;

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'spoilers', spoilers, 'plot_flag', plot_flag);


%%  %%%% Example: 2 site pyr->lac, pyruvate: 2D GRE, lactate: 2D GRE  %%%%%%
clear; close all;
%set parameters

%flip angles
flips.P = 20;
flips.L = 30;

% TR
TRP = .125; TRL = .125; %s 
TR.P = TRP * 16;
TR.L = TRL * 16;

%relaxation
R1 = [1/30, 1/25];
k = [0.002];

TempRes = 4;
Nt= 30;
Mz0 = [1 0.1];
spoilers = [1 1];
plot_flag = 1;

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'spoilers', spoilers, 'plot_flag', plot_flag);


%% %%%% Example: 2 site pyr->lac W/INPUT, pyruvate: 2D GRE, lactate: 3D bSSFP  %%%%%%

%set parameters
clear; close all;

%flip angles
flips.P = 20;
FAL = 60;
Nz = 16; interleaves = 4;
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60];

% TR
TRP = .18; TRL = .01529; %s 
TR.P = TRP * 16;
TR.L = repmat(TRL, [1, Nz*interleaves]);
cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL];

%relaxation
R1 = [1/30, 1/25];
R2.L = 1/1;
k = [0.005]; %[0.02];

TempRes = 4;
Nt= 30;
Mz0 = [0 0];
spoilers = [1 0];
plot_flag = 1;
scales = [1 6.905];

%input function
Tarrival = 0; Tbolus = 8;
[input_function, ~] = realistic_input_function(Nt, TempRes, Tarrival, Tbolus);

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'spoilers', spoilers, 'scales', scales, 'input', input_function, 'plot_flag', plot_flag);



%% %%%% Example: 2 site pyr->lac W/ INPUT, pyruvate: 3D GRE, lactate: 3D bSSFP  %%%%%%

%set parameters
clear; close all;

%flip angles
FAP = 3;
FAL = 60;
Nz = 16; interleaves = 4;
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60];

% TR
TRP = .18; TRL = .01529; %s 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz*interleaves]);
%cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL];
cat_TR = [TRL/2, TRL/2, TRL/2, TRL/2, TRL/2, TRL/2];

%relaxation
R1 = [1/30, 1/25];
R2.L = 1/1;
k = [0.002];

TempRes = 4;
Nt= 30;
Mz0 = [0 0];
spoilers = [1 0];
scales = [1 1];
plot_flag = 1;

%input function
Tarrival = 0; Tbolus = 8;
[input_function, ~] = realistic_input_function(Nt, TempRes, Tarrival, Tbolus);

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'spoilers', spoilers, 'scales', scales, 'input', input_function, 'plot_flag', plot_flag);

%% %%%% Example: 2 site pyr->lac W/ INPUT, pyruvate: 3D GRE, lactate: 3D GRE  %%%%%%

%set parameters
clear; close all;

%flip angles
FAP = 3;
FAL = 7.67;
Nz = 16;
flips.P = repmat(FAP, [1, Nz]);
flips.L = repmat(FAL, [1, Nz]);

% TR
TRP = .125; TRL = .125; %s 
TR.P = repmat(TRP, [1, Nz]);
TR.L = repmat(TRL, [1, Nz]);

%relaxation
R1 = [1/30, 1/25];
k = [0.002]; 

TempRes = 4;
Nt= 30;
Mz0 = [0 0];
spoilers = [1 1];
plot_flag = 1;
Rinj = 0.1;

%input function
Tarrival = 0; Tbolus = 8;
[input_function, ~] = realistic_input_function(Nt, TempRes, Tarrival, Tbolus);

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'spoilers', spoilers, 'input', input_function, 'plot_flag', plot_flag);



%% %%%% Example: 4 site pyr->lac&bic&ala W/ INPUT, pyruvate: 2D GRE, lactate: 3D bSSFP, bicarb: 2D GRE, alanine: 3D bSSFP  %%%%%%

%set parameters
clear; close all;

TempRes = 5;
Nt= 30;

%flip angles
flips.P = 20;
FAL = 60;
Nz = 16; interleaves = 4;
flips.L = repmat(FAL, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);
cat_flips = [4, -16, 24, -36, 48, -60];
flips.B = 30;
FAA = 60;
flips.A = repmat(FAA, [1, Nz*interleaves]) .* repmat([1, -1], [1, Nz*interleaves/2]);

% TR
TRP = .08; TRL = .01529; TRB = .08; TRA = .01529; %s 
TR.P = TRP * 16;
TR.L = repmat(TRL, [1, Nz*interleaves]);
TR.B = TRB * 16;
TR.A =  repmat(TRA, [1, Nz*interleaves]);

cat_TR = [TRL, TRL, TRL, TRL, TRL, TRL];

%relaxation
R1 = [1/30, 1/25, 1/25, 1/25];
R2.L = 1/1; R2.A = 1/1;
k = [0.02 0.01 0.005];

Mz0 = [0 0 0 0];
spoilers = [1 0 1 0];
scales = [1 1 1 1];
plot_flag = 1;

%input function
Tarrival = 0; Tbolus = 8;
[input_function, ~] = realistic_input_function(Nt, TempRes, Tarrival, Tbolus);

%run model
[Mxy, Mz] = sim_multisite_bSSFP(flips, TR, TempRes, R1, k, Nt, Mz0, 'R2', R2, 'cat_flips', cat_flips, 'cat_TR', cat_TR, 'spoilers', spoilers, 'scales', scales, 'input', input_function, 'plot_flag', plot_flag);

