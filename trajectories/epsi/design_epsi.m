function [g, ktraj_g, gparams, opts] = design_epsi(epsi_type, ramp_sampling, spatial_res, spatial_fov, spec_res, spec_bw, opts)
% DESIGN_EPSI - design an echo-planar spectroscopic imaging gradient
%
% [g, ktraj_g, gparams] = design_epsi(epsi_type, ramp_sampling, spatial_res, spatial_fov, spec_res, spec_bw, [opts])
%
% INPUTS
%	epsi_type - 'flyback' (default), 'symmetric'
%	ramp_sampling - 1=on, 0=off (partial ramp sampling also supported by
%       values between 0-1)
%   spatial_res - Spatial resolution (cm)
%   spatial_fov - Spatial FOV (cm), maybe increased by design
%   spec_res - Spectral resolution (Hz), assumes half-echo (half-Fourier) time sampling
%   spec_bw - Spectral bandwidth (Hz).  maybe increased design
%	opts (optional) - structure defining options for the EPSI design.  This includes:
%       'max_slew' (default = 20 G/cm/ms), 'max_g' (default = 5 G/cm/ms),
%       'samp_rate' - waveform sampling rate (default = 4e-6 s),
%       'GAMMA' (13C is default = 1071 Hz/G)
% OUTPUTS
%   g - gradient (G/cm)
%   ktraj_g - k-space trajectory (1/cm)
%   gparams - structure containing other EPSI parameters:
%       'data_samp_rate' - data sampling rate (s)
%       'sampling_delay' - delay between first gradient point and first
%           readout sample
%       'n_read' - samples to readout per gradient lobe, at data_samp_rate
%       'n_skip' - samples to skip between readout sections, at data_samp_rate
%       'Nlobes' - number of readout gradient lobes
%       'spatial_fov' - actual spatial FOV (cm)
%       'spatial_res' - actual spatial resolution (cm)
%       'spec_res' - actual spectral resolution (Hz)
%       'spec_bw' - actual spectral bandwidth (Hz)
%       'epsi_type' - 'flyback', 'symmetric'
%       'ramp_sampling' - fraction of the ramps used for recon
%       'samp_rate' - waveform sampling rate
%       'GAMMA' - gyromagnetic ratio


%
% Author: Peder E. Z. Larson
%
% (c)2014-2017 The Regents of the University of California.
% All Rights Reserved.

% Note that 'symmetric' is Nyquist-sampled symmetric, meaning nowhere in k-space is
% the Nyquist criteria violated, and there will also be oversampling.


if nargin < 7
    opts = [];
end
if ~isfield(opts, 'GAMMA')
    opts.GAMMA = 1071;  % Hz/G
end
if ~isfield(opts, 'max_g')
    opts.max_g = 5;  % G/cm
end
if ~isfield(opts, 'max_slew')
    opts.max_slew = 20; % G/cm/ms
end
if ~isfield(opts, 'samp_rate')
    opts.samp_rate = 4e-6; % s
end


% design gradient lobes
% dzg will determine whether spatial_res and spec_bw specifications are
% feasible
switch epsi_type
    case 'symmetric'
        [g_readlobe, n_samp_delay, gparams.data_samp_rate, gparams.n_skip, gparams.n_read, n_plateau_out, n_ramp_out] = ...
            dzg_read(floor(1/(2*spec_bw)/opts.samp_rate)*opts.samp_rate, 1/spatial_res, ramp_sampling);
        g_dephase = dzg_short(sum(g_readlobe)/2 * opts.samp_rate * opts.GAMMA);
        Nlobes = ceil(1/(2*spec_res * 2*length(g_readlobe)* opts.samp_rate));
        g = [-g_dephase, kron(ones(1,Nlobes),[g_readlobe -g_readlobe])];
        
        gparams.pw_read_plateau = n_plateau_out*opts.samp_rate;
        gparams.pw_read_ramp = n_ramp_out*opts.samp_rate;
        gparams.Nlobes = 2*Nlobes;
        gparams.spec_res = 1/(2*Nlobes*2*length(g_readlobe)* opts.samp_rate);
        gparams.spec_bw = 1 / (2*length(g_readlobe)*opts.samp_rate);
        
    case 'flyback'
        % prioritizes spec_bw
        %  g_fblobe = dzg_short(1/spatial_res);
        %  [g_readlobe gparams.data_samp_rate, gparams.n_skip, gparams.n_read] = ...
        %       dzg_read(1/spec_bw - length(g_fblobe)*opts.samp_rate, 1/spatial_res, ramp_sampling);
        
        % prioritizes spatial_res
        [g_readlobe, g_fblobe, n_samp_delay, gparams.data_samp_rate, gparams.n_skip, gparams.n_read, n_plateau1_out, n_plateau2_out, n_ramp1_out, n_ramp2_out] = ...
            dzg_flyback(floor(1/spec_bw/opts.samp_rate)*opts.samp_rate, 1/spatial_res, ramp_sampling,opts, spatial_fov);
        
        g_dephase = dzg_short(sum(g_readlobe)/2 * opts.samp_rate * opts.GAMMA);
        Nlobes = ceil(1/(2*spec_res * (length(g_readlobe)+length(g_fblobe))* opts.samp_rate));
        g = [-g_dephase, kron(ones(1,Nlobes),[g_readlobe -g_fblobe])];
        
        gparams.pw_read_plateau = n_plateau1_out*opts.samp_rate;
        gparams.pw_read_ramp = n_ramp1_out*opts.samp_rate;
        gparams.pw_fb_plateau = n_plateau2_out*opts.samp_rate;
        gparams.pw_fb_ramp = n_ramp2_out*opts.samp_rate;
        gparams.Nlobes = Nlobes;
        gparams.spec_res = 1/(2*Nlobes*(length(g_readlobe)+length(g_fblobe))* opts.samp_rate);
        gparams.spec_bw = 1 / ((length(g_readlobe)+length(g_fblobe))*opts.samp_rate);
        
        
end


gparams.sampling_delay = (length(g_dephase) + n_samp_delay) * opts.samp_rate;
gparams.spatial_fov = 1/ (opts.GAMMA * max(g_readlobe) * gparams.data_samp_rate);
gparams.spatial_res = spatial_res;
gparams.epsi_type = epsi_type;
gparams.ramp_sampling = ramp_sampling;
gparams.samp_rate = opts.samp_rate;
gparams.GAMMA = opts.GAMMA;
% Add SNR efficiency estimate?

ktraj_g = opts.GAMMA*opts.samp_rate*cumsum(g);

% end of main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nested functions:

    function [g, n_samp_delay, data_samp_rate, n_skip, n_read, n_plateau, n_ramp] = dzg_read(T, Ak, framp)
        % T - gradient duration (s)
        % Ak - total gradient area for data sampling (1/cm)
        % framp - fraction of ramp samples to use, [0,1]
        
        S = opts.max_slew*1e3;
        A = Ak / opts.GAMMA;
        
        t_ramp = ( T*S - sqrt( (T*S)^2 - 4*(2-framp)*S*A) )/(2*(2-framp)*S);
        if ~isreal(t_ramp)
            error('Cannot design EPSI for chosen spatial_res, spec_bw & ramp_sampling.  Try relaxing these contraints.')
        end
        g_ideal = S*t_ramp;
        
        % data sampling rate
        data_samp_rate_factor = floor(1 / (opts.GAMMA * g_ideal * spatial_fov) / opts.samp_rate);
        data_samp_rate = data_samp_rate_factor * opts.samp_rate;
        
        % only stretch gradient to maintain within slew, g limits
        % in units of opts.samp_rate
        if framp == 0
            % enforce data_samp_rate on all gradient sections for
            % simplifying data reconstruction
            n_ramp = ceil(t_ramp / data_samp_rate ) * data_samp_rate_factor ;
            n_plateau = floor( (T - 2*n_ramp*opts.samp_rate)/data_samp_rate ) * data_samp_rate_factor ;
        else
            % only enforce data_samp_rate over total gradient length
            n_ramp = ceil(t_ramp / opts.samp_rate );
            n_plateau = floor( T/data_samp_rate)* data_samp_rate_factor  - 2*n_ramp ;
            
        end
        
        n_samp_delay = floor(n_ramp * (1-framp));
        
        % in units of data_samp_rate
        n_read = (n_plateau+2*(n_ramp-n_samp_delay)) / data_samp_rate_factor;
        n_skip = 2*n_samp_delay / data_samp_rate_factor;
        
        g_max = A / ( (n_plateau+(n_ramp-n_samp_delay))*opts.samp_rate);
        
        g = [ [0.5:n_ramp-.5]/n_ramp ones(1, n_plateau) [n_ramp-0.5:-1:0.5]/n_ramp ] * g_max;
        
        while max(abs(diff(g))/opts.samp_rate) > S
            % error('Slew rate exceeded')
            % if slew rate violated, increase ramp time, but keep total
            % duration the same (not sure this will work in all conditions
            % though)
            if framp == 0
                n_ramp = n_ramp + data_samp_rate_factor ;
                n_plateau = floor( (T - 2*n_ramp*opts.samp_rate)/data_samp_rate ) * data_samp_rate_factor ;
                
            else
                n_ramp = n_ramp+1;
                n_plateau = floor( T/data_samp_rate)* data_samp_rate_factor  - 2*n_ramp ;
                
            end
            n_samp_delay = floor(n_ramp * (1-framp));
            
            % in units of data_samp_rate
            n_read = (n_plateau+2*(n_ramp-n_samp_delay)) / data_samp_rate_factor;
            n_skip = 2*n_samp_delay / data_samp_rate_factor;
            
            g_max = A / ( (n_plateau+(n_ramp-n_samp_delay))*opts.samp_rate);
            % need to recheck slew here, if g_max> g_ideal?
            
            
            g = [ [0.5:n_ramp-.5]/n_ramp ones(1, n_plateau) [n_ramp-0.5:-1:0.5]/n_ramp ] * g_max;
        end
        
        if g_max > opts.max_g
            error('Cannot design EPSI for chosen spatial_res, spec_bw & ramp_sampling.  Try relaxing these contraints.')
        end
        
        
        
    end

    function g = dzg_short(Ak)
        % Ak - total gradient area (1/cm)
        S = opts.max_slew*1e3;
        A = Ak / opts.GAMMA;
        
        gmax_ideal = sqrt( A * S ); % without gradient limitation
        n_ramp = ceil(gmax_ideal / (S * opts.samp_rate));
        n_ramp_max = ceil(opts.max_g / (S * opts.samp_rate));
        if n_ramp > n_ramp_max
            n_ramp = n_ramp_max;
            n_plateau = ceil( (A - opts.max_g* n_ramp*opts.samp_rate) / (opts.max_g * opts.samp_rate));
        else
            n_plateau = 0;
        end
        
        g_max = A / ( (n_plateau+n_ramp)*opts.samp_rate);
        g = [ [0.5:n_ramp-.5]/n_ramp ones(1, n_plateau) [n_ramp-0.5:-1:0.5]/n_ramp ] * g_max;
        
        if max(abs(diff(g))/opts.samp_rate) > S
            error('Slew rate exceeded')
        end
    end
end

function [gread, gfb, n_samp_delay, data_samp_rate, n_skip, n_read, n_plateau1, n_plateau2, n_ramp1, n_ramp2] = dzg_flyback(T, Ak, framp,opts,spatial_fov)
% T - total gradient duration (s)
% Ak - total gradient area for data sampling (1/cm)
% framp - fraction of ramp samples to use, [0,1]

S = opts.max_slew*1e3;
A = Ak / opts.GAMMA;
n_ramp_max = ceil(opts.max_g / (S * opts.samp_rate));

syms d1 d2 D1 D2 positive

% determine whether plateau samples required on flyback lobe
solns = solve(A == (D1 + framp*d1)*S*d1, T == D1 + 2*(d1+d2), S*d1*(D1+d1) == S*d2^2, 'Real', true);
if isempty(solns)
    error('Cannot design EPSI for chosen spatial_res, spec_bw & ramp_sampling.  Try relaxing these contraints.')
end
SD1 = eval(solns.D1); SD2 = zeros(size(SD1)); Sd1 = eval(solns.d1); Sd2 = eval(solns.d2);  % This may require Simulink toolbox, use 'float' instead?

n_ramp1 = ceil(Sd1 / opts.samp_rate);
n_ramp2 = ceil(Sd2 / opts.samp_rate);
n_plateau1 = floor( (T - 2*(n_ramp1+n_ramp2)*opts.samp_rate)/opts.samp_rate );
n_plateau2 = zeros(size(n_plateau1));

Isoln = find(n_ramp1 <= n_ramp_max & n_ramp2 <= n_ramp_max);
if isempty(Isoln)
    % plateau samples required
    solns = solve(A == (D1 + framp*d1)*S*d1, T == D1 + D2 + 2*(d1+d2), S*d1*(D1+d1) == S*d2*(D2+d2), opts.max_g == S*d2, 'Real', true);
    if isempty(solns)
        error('Cannot design EPSI for chosen spatial_res, spec_bw & ramp_sampling.  Try relaxing these contraints.')
    end
    SD1 = eval(solns.D1); SD2 = eval(solns.D2); Sd1 = eval(solns.d1); Sd2 = eval(solns.d2);
    n_ramp1 = ceil(Sd1 / opts.samp_rate);
    n_ramp2 = ceil(Sd2 / opts.samp_rate);
    n_plateau1 = floor( (T - SD2 - 2*(n_ramp1+n_ramp2)*opts.samp_rate)/opts.samp_rate );
    n_plateau2 = ceil( (T - (n_plateau1 + 2*(n_ramp1+n_ramp2))*opts.samp_rate)/opts.samp_rate );
    Isoln = find(n_ramp1 <= n_ramp_max & n_ramp2 <= n_ramp_max);
end

if length(Isoln > 1)
    % Choose solution with longest readout plateau duration
    [temp ImaxD1] = max(SD1(Isoln));
    Isoln = Isoln(ImaxD1);
end
n_ramp1 = n_ramp1(Isoln); n_ramp2 = n_ramp2(Isoln);
n_plateau1 = n_plateau1(Isoln);    n_plateau2 = n_plateau2(Isoln);

% correct waveform to use a specific data sampling rate
g_ideal = S*Sd1(Isoln);
data_samp_rate_factor = floor(1 / (opts.GAMMA * g_ideal * spatial_fov) / opts.samp_rate);
data_samp_rate = data_samp_rate_factor * opts.samp_rate;

% only stretch gradient to maintain within slew, g limits
% in units of opts.samp_rate
if framp == 0
    % enforce data_samp_rate on read gradient ramps and flyback lobe for
    % simplifying data reconstruction
    n_ramp1 = ceil(Sd1(Isoln) / data_samp_rate ) * data_samp_rate_factor;
    n_ramp2 = ceil(Sd2(Isoln) / opts.samp_rate );
    n_plateau1 = floor( (T - SD2(Isoln) - 2*(n_ramp1+n_ramp2)*opts.samp_rate)/data_samp_rate ) * data_samp_rate_factor;
    n_plateau2 = ceil( (T - (n_plateau1 + 2*n_ramp1)*opts.samp_rate)/data_samp_rate)* data_samp_rate_factor  - 2*n_ramp2;
else
    % enforce data_samp_rate over read and flyback gradient lobes
    n_plateau1 = floor( (SD1(Isoln) + 2*Sd1(Isoln))/data_samp_rate)* data_samp_rate_factor  - 2*n_ramp1;
    n_plateau2 = ceil( (SD2(Isoln) + 2*Sd2(Isoln))/data_samp_rate)* data_samp_rate_factor  - 2*n_ramp2;
end

n_samp_delay = floor(n_ramp1 * (1-framp));

% in units of data_samp_rate
n_read = (n_plateau1+2*(n_ramp1-n_samp_delay)) / data_samp_rate_factor;
n_skip = (2*n_samp_delay + n_plateau2 + 2*n_ramp2) / data_samp_rate_factor;

g_max1 = A / ( (n_plateau1+(n_ramp1-n_samp_delay))*opts.samp_rate);

g_max2 = g_max1 * (n_plateau1 + n_ramp1) / (n_plateau2 + n_ramp2);

if g_max1 > opts.max_g || g_max2 > opts.max_g
    error('Cannot design EPSI for chosen spatial_res, spec_bw & ramp_sampling.  Try relaxing these contraints.')
end

gread = [ [0.5:n_ramp1-.5]/n_ramp1 ones(1, n_plateau1) [n_ramp1-0.5:-1:0.5]/n_ramp1 ] * g_max1;
gfb = [ [0.5:n_ramp2-.5]/n_ramp2 ones(1, n_plateau2) [n_ramp2-0.5:-1:0.5]/n_ramp2 ] * g_max2;

if max(abs(diff([gread,-gfb]))/opts.samp_rate) > S
    error('Slew rate exceeded')
end

end



