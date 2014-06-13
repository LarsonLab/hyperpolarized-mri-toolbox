% EPSI design examples
%
% Author: Peder E. Z. Larson
%
% (c)2014 The Regents of the University of California.
% All Rights Reserved.

validate = 0;

%% Design for a clinical 3T MRI

% Note that default gradient options are for 13C on a typical clinical MRI system

spatial_res = 0.3; % cm
spatial_fov = 21.9; % cm
spec_res = 5; % Hz
spec_bw = 470; % Hz - cover from pyruvate to lactate

epsi_type = 'symmetric'; ramp_sampling = 1;

[g, ktraj_g, gparams, opts] = design_epsi(epsi_type, ramp_sampling, spatial_res, spatial_fov, spec_res, spec_bw);

disp(gparams)

figure
subplot(211), plot([1:length(g)]*opts.samp_rate, g), ylabel('Gradient (G/cm)')
subplot(212), plot([1:length(g)]*opts.samp_rate, ktraj_g), ylabel('k-space (1/cm)'), xlabel('time (s)')

if validate
    gnew = g; ktraj_gnew = ktraj_g;
    load test_data/epsi_3t_example
    if max(abs(gnew - g)) > eps || max(abs(ktraj_gnew - ktraj_g)) > eps
        error('Validation failed for 3T EPSI example')
    else
        disp('Validation passed for 3T EPSI example')
    end
end
pause

%% Design for a high performance, small-animal 1T MRI system
clear opts
opts.max_g = 45; % G/cm
opts.max_slew = 45 / .25; % G/cm/ms
opts.samp_rate = 2e-6; % s  (system allows for 0.1 us, but this is much more than necessary)
opts.GAMMA = 1071; % Hz/G

spatial_res = 0.2; % cm
spatial_fov = 16; % cm
spec_res = 5; % Hz
spec_bw = 275; % Hz - should cover from bicarbonate to lactate

epsi_type = 'flyback'; ramp_sampling = 0;  % high gradient strengths should allow for use of flyback without ramp sampling

[g, ktraj_g, gparams] = design_epsi(epsi_type, ramp_sampling, spatial_res, spatial_fov, spec_res, spec_bw, opts);

disp(gparams)

figure
subplot(211), plot([1:length(g)]*opts.samp_rate, g), ylabel('Gradient (G/cm)')
subplot(212), plot([1:length(g)]*opts.samp_rate, ktraj_g), ylabel('k-space (1/cm)'), xlabel('time (s)')

if validate
    gnew = g; ktraj_gnew = ktraj_g;
    load test_data/epsi_1t_example
    if max(abs(gnew - g)) > eps || max(abs(ktraj_gnew - ktraj_g)) > eps
        error('Validation failed for 1T EPSI example')
    else
        disp('Validation passed for 1T EPSI example')
    end
end
