function write_epsi(g, gparams, filename, format)
% WRITE_EPSI - write an echo-planar spectroscopic imaging gradient
%   companion to design_epsi() function.
%
% write_epsi(
%
% INPUTS
%   g - gradient (G/cm)
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
%
%
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

%
% Author: Peder E. Z. Larson
%
% (c) 2017 The Regents of the University of California.
% All Rights Reserved.

if (nargin < 3) || isempty(filename)

filename = input('Root file name: (leave empty to not save) ', 's');
if isempty(filename)
    fprintf(1,'Not saving files \n');
    return;
end;
end

if (nargin < 4) || isempty(format)
    format = 'GE';
else
    switch format,
        case {'GE', 'Varian'}
        otherwise
            error('Format save type of: %s not recognized', format);
    end;
end

dat_name = sprintf('%s.dat', filename);
fid = fopen(dat_name, 'w');
if fid == -1,
    fprintf(1, 'Error opening %s \n', dat_name);
    return;
end;

fprintf(fid, 'EPSI waveform %s, flyback = %d\n', filename, strcmp(gparams.epsi_type,'flyback'));
fprintf(fid, '%d #pw(us)\n', length(g)*gparams.samp_rate*1e6);
fprintf(fid, '%d #res\n', length(g));
fprintf(fid, '%f #spectral_width(Hz)\n', gparams.spec_bw);
fprintf(fid, '%d #filter_bandwidth(Hz)\n', round(1/gparams.data_samp_rate));

fprintf(fid, '%6d #filter_points\n', ceil(gparams.Nlobes * (gparams.n_read + gparams.n_skip) /2)*2); %must be even

fprintf(fid, '%f #nom_gradient(G/cm)\n', -max(abs(g)));

fprintf(fid, '%f #nom_res(mm)\n', gparams.spatial_res*10);
fprintf(fid, '%d #sampling_delay(us)\n', gparams.sampling_delay*1e6);
fprintf(fid, '%d #nlobes\n', gparams.Nlobes);

fprintf(fid, '%d #offset_pts(ramp_points)\n', round(gparams.pw_read_ramp/gparams.samp_rate)); %??

fprintf(fid, '%d #nspatial\n', floor(gparams.spatial_fov / gparams.spatial_res));
fprintf(fid, '%d #nskip\n', gparams.n_skip);
% change below to integers?
fprintf(fid, '%f #dwell_time(us)\n', gparams.data_samp_rate*1e6);
fprintf(fid, '%f #pw_read_ramp(us)\n', gparams.pw_read_ramp*1e6);
fprintf(fid, '%f #pw_read_plateau(us)\n', gparams.pw_read_plateau*1e6);
if strcmp(gparams.epsi_type, 'flyback')
    fprintf(fid, '%f #pw_fb_ramp(us)\n', gparams.pw_fb_ramp*1e6);
    fprintf(fid, '%f #pw_fb_plateau(us)\n', gparams.pw_fb_plateau*1e6);
else
    fprintf(fid, '%f #pw_read_ramp(us)\n', gparams.pw_read_ramp*1e6); % for positive and negative lobe durations
    fprintf(fid, '%f #pw_read_plateau(us)\n', gparams.pw_read_plateau*1e6);
end


fclose(fid);



switch (format)
    case 'GE'
        
        % write out gradient file
        signa(g,[filename '.wav'],1);
        
    case 'Varian'
        
        fid = fopen(sprintf('%s.GRD', filename),'wt');
        fprintf(fid,'# %s\n', sprintf('%s.GRD', filename));
        fprintf(fid,'# ***************************************************\n');
        fprintf(fid,'# EPSI Waveform\n');
        fprintf(fid,'# see %s.dat for waveform parameters\n',filename);
        fprintf(fid,'# ***************************************************\n');
        fprintf(fid,'%d  1\n',round(32767*g/max(abs(g))));
        fclose(fid);
        
        
end

end
