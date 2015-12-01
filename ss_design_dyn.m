function [g, rf, fs_best, z_plot, f_plot, m_plot] = ...
    ss_design_dyn(z_thk, z_tb, z_de, f, a_angs, dr,...
    ptype, z_ftype, s_ftype, ss_type, ...
    f_off, dbg, no_plot);


% SS_DESIGN_DYN - Design spectral-spatial pulse
%   Fixing some parameters for dynamic spectral-spatial designs
%
% [g, rf, fs_best, z_plot, f_plot, m_plot] = ...
%     ss_design(z_thk, z_tb, z_de, f, a_angs, de,...
%               ptype, z_ftype, s_ftype, ss_type, ...
%               f_off, dbg, no_plot)
%
% INPUTS
%   z_thk - slice thickness (cm)
%   z_tb - spatial time-bandwidth
%   z_de - spatial ripples, [pass_ripple, stop_ripple]
%   f - spectral band edge specification (Hz)
%   a_angs - spectral band flip angle specification (radians)
%   de - relative spectral band ripple
%   ptype - spatial pulse type: 'ex' (default), 'se', 'sat', 'inv'
%   z_ftype - spatial filter type: 'ms', 'ls', 'pm' (default), 'min', 'max'
%   s_ftype - spectral filter type: 'min' (default), 'max', 'lin'
%   ss_type - spectral-spatial type: 'Flyback Whole' (default),
%       'Flyback Half', 'EP Whole', 'EP Half',
%       'EP Whole Opp-Null', 'EP Half Opp-Null'
%   f_off - center frequency (empty to let ss_design choose)
%   dbg - print debug messages: 0-none (default), 1-some, 2-more
%   no_plot - set to 1 to turn off plotting
%
% OUTPUTS
%   g - gradient (G/cm)
%   rf - RF (G)
%   fs_best - spectral sampling frequency (Hz)
%   z_plot, f_plot, m_plot - ss_plot() outputs, ranges and values in figure
%
% See scripts in examples/ folder demonstrations of how to use this
% function.

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

%
% $Header: /home/adam/cvsroot/src/ss/ss_design.m,v 1.33 2014/05/22 22:44:36 peder Exp $
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check all inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (nargin < 6),
    error(['Usage: ss_design(z_thk, z_tb, z_d, f, a_angs, d,' ...
        'ptype, z_fttype, s_ftype, ss_type, foff)']);
end;

% Check ptype
%
if (nargin < 7) || isempty(ptype),
    ptype = 'ex';
else
    switch ptype,
        case {'ex', 'se', 'sat', 'inv'}
        otherwise
            error(sprintf(['Spatial pulse type (ptype) of: %s not' ...
                ' recognized'], ptype));
    end;
end;

% Check z_ftype
%
if (nargin < 8) || isempty(z_ftype),
    z_ftype = 'pm';
else
    switch z_ftype,
        case {'ms', 'ls', 'pm', 'min', 'max'}
        otherwise
            error(sprintf(['Spatial filter type (z_ftype) of: %s not' ...
                ' recognized'], z_ftype));
    end;
end;

% Check s_ftype
%
if (nargin < 9) || isempty(s_ftype),
    s_ftype = 'min';
else
    switch s_ftype,
        case {'min', 'max', 'lin'}
        otherwise
            error(sprintf(['Spectral filter type (s_ftype) of: %s not' ...
                ' recognized'], s_ftype));
    end;
end;

% Check ss_type
%
if (nargin < 10) || isempty(ss_type),
    ss_type = 'Flyback Whole';
else
    switch ss_type,
        case {'Flyback Whole', 'Flyback Half', ...
                'EP Whole', 'EP Half', 'EP Whole Opp-Null', ...
                'EP Half Opp-Null'}
        otherwise
            error(sprintf(['Spectral-spatial type (ss_type) of: %s not' ...
                ' recognized'], ss_type));
    end;
end;

% Check f_off
%
if (nargin < 11) || isempty(f_off),
    f_off = [];
end;

% Check dbg
%
if (nargin < 12) || isempty(dbg),
    dbg = 0;
end;

if (nargin < 13) || isempty(no_plot),
    no_plot = 0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize globals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss_globals;

global SS_MAX_DURATION SS_MIN_ORDER;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert "a_angs" to Beta polynomial "a"
% Convert effective ripples to polynomial ripples
% depending on pulse type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = size(a_angs,2);

a = zeros(size(a_angs));
for t = 1:Nt
    ang(t) = max(a_angs(:,t));
    
    if SS_SLR_FLAG,
        a(:,t) = sin(a_angs(:,t)/2)/sin(ang(t)/2);
        %a = asin(sin(ang)*a)/ang;
    else
        a(:,t) = a_angs(:,t)/ang(t);
    end
    de(:,t) = a_angs(:,t) *dr;
    
    [d(:,t),a(:,t),ang(t)] = rf_ripple(de(:,t), a(:,t), ang(t), ptype);
    z_d(:,t) = rf_ripple(z_de, [1 0], ang(t), ptype);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate gradient lobe required to achieve spatial
% prescription
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate cycles/cm required
%
kz_max = z_tb / z_thk;		% cycles/cm
kz_area = kz_max / SS_GAMMA;	% G/cm * s

% Note:  If ss_type is EP, then SS_EQUAL_LOBES must be 1
%
if strfind(ss_type, 'EP')
    SS_EQUAL_LOBES = 1;
else;
    SS_EQUAL_LOBES = 0;
end;

% Get SS gradient lobes that gives this area with
% the desired fraction of sloped gradient included
% (this is to limit amount of versing required)
%
[gpos, gneg, g1, g2, g3] = ...
    grad_ss(kz_area, [], SS_VERSE_FRAC, SS_MXG, SS_MXS, ...
    SS_TS, SS_EQUAL_LOBES);

% Calculate maximum spectral sampling frequency
%
lobe = [gpos gneg];
t_poslobe = length(gpos) * SS_TS;
t_lobe = length(lobe) * SS_TS;
fs_scale = 1;
if strfind(ss_type, 'EP')
    fs_max = 2/t_lobe;
    if ~strfind(ss_type, 'Opp-Null')
        fs_scale = 1/2;		% If true null, make sure
    end;				% aliasing check is for half-frequency
else
    fs_max = 1/t_lobe;
end;

% Estimate sampling frequencies that won't cause
% incompatible aliasing of stopbands/passbands
%
switch ss_type,
    case {'Flyback Half', 'EP Half', 'EP Half Opp-Null'},
        sym_flag = 1;
    otherwise
        sym_flag = 0;
end;

fdiff = diff(f);
fwidths = sort(fdiff(1:2:end), 2, 'descend');
fs_min = sum(fwidths(1:2))/2;		% Loose estimate on lower bound
df = (fs_max-fs_min)/(SS_NUM_FS_TEST-1);

Nbands = length(f)/2;
a_test = 1:Nbands;
% test to see if flip angles are the same for all time, in which case they
% can alias on top of each other
for I1 = 1:Nbands-1
    for I2 = I1+1:Nbands
        if all(a_angs(I1,:) == a_angs(I2,:))
            a_test(I2) = a_test(I1);
        end
    end
end
d_test = .01*ones(1,Nbands);

for idx = 1:SS_NUM_FS_TEST,
    fs_test(idx) = fs_min + (idx-1)*df;
    % Make sure fs_cur is integer number of samples
    %
    nsamp = ceil(1/(fs_test(idx)*SS_TS));
    fs_test(idx) = 1/(nsamp*SS_TS);
    
    % Test aliasing of frequencies into operating BW
    %
    [f_a, a_a, d_a, fo] = ss_alias(f,a_test,d_test,f_off,fs_test(idx)*fs_scale,sym_flag);
    if (isempty(f_a)),
        fs_ok(idx) = 0;
    else
        fs_ok(idx) = 1;
    end;
end;
if (dbg >= 2)
    figure;
    plot(fs_test, fs_ok, '*');
    title('Feasible Spectral Sampling Frequencies');
    xlabel('Frequency [Hz]');
    ylabel('Flag (1=OK, 0=No Good)');
    drawnow;
    pause(1);
end;

if (~any(fs_ok == 1)),
    fs_over = fs_test(end);
    while isempty(f_a)
        fs_over = fs_over + df;
        [f_a, a_a, d_a, fo] = ss_alias(f,a_test,d_test,f_off,fs_over*fs_scale,sym_flag);
    end;
    
    fprintf(1,'ss_design: Incompatible aliasing of frequency spec\n');
    fprintf(1,'at all tested frequencies. \n');
    fprintf(1,'Current max sampling is: %6.1f\n', fs_max);
    fprintf(1,'Estimated required sampling is: %6.1f\n', fs_over);
    fprintf(1, 'Try any of the following to incr:\n');
    fprintf(1,'   - Decrease spatial TBW\n');
    fprintf(1,'   - Increase slice thickness\n');
    fprintf(1,'   - Increase VERSE fraction\n');
    fprintf(1,'   - Decrease frequency band widths\n');
    switch (ss_type),
        case {'Flyback Whole', 'Flyback Half'}
            fprintf(1,'   - Try EPI type SS\n');
    end;
    error('No good fs');
end;
fs_bands_left = find(diff([0 fs_ok]) == 1);
fs_bands_right = find(diff([fs_ok 0]) == -1);
fs_bands = [fs_bands_left; fs_bands_right];
fs_bands = fs_bands(:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate on lobe width trying to meet B1 requirements
% with minimum-time pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nsolutions = 0;

if (SS_NUM_LOBE_ITERS == 1),
    fs_top = fs_bands(end);
    fs_best = fs_test(fs_top);
    
    % Call design script for SS that calculates rf, g that meets
    % sampling requirements
    %
    b1_set = 0;  dur_set = 0; pow_set = 0;         nsolutions = 1;
    pw_set = 0;
    
    % First, find what sampling frequencies work and required pulse
    % durations
    for t = 1:Nt
        
        if strfind(ss_type, 'Flyback')
            [rftemp{t}, gtemp{t}] = ss_flyback(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs_best, ...
                ptype, z_ftype, s_ftype, ss_type, ...
                f_off, dbg);
        else
            [rftemp{t}, gtemp{t}] = ss_ep(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs_best, ...
                ptype, z_ftype, s_ftype, ss_type, ...
                f_off, dbg);
        end;
        
        if ~isempty(rftemp{t}),
            pw_set = max([length(rftemp{t}), pw_set]);
        else
            nsolutions = 0;
        end;
    end
    
    if nsolutions == 1
        
        % fix pulse duration for all RF pulses
        SS_MAX_DURATION = pw_set*SS_TS;  SS_MIN_ORDER = 0;
        for t = 1:Nt
            
            if strfind(ss_type, 'Flyback')
                [rftemp{t}, gtemp{t}] = ss_flyback(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs_best, ...
                    ptype, z_ftype, s_ftype, ss_type, ...
                    f_off, dbg);
            else
                [rftemp{t}, gtemp{t}] = ss_ep(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs_best, ...
                    ptype, z_ftype, s_ftype, ss_type, ...
                    f_off, dbg);
            end;
            
            if ~isempty(rftemp{t}),
                b1_set = max([abs(rftemp{t}),b1_set]);
                dur_set = max([length(rftemp{t}) * SS_TS, dur_set]);
                pow_set = max([sum(abs(rftemp{t}).^2) * SS_TS, pow_set]);
            else
                nsolutions = 0;
            end;
            
        end
        
    end
    
    if nsolutions > 0
        fprintf(1,'Solution(s) exists!\n');
        
        % fix rephasing gradient area to be same as initial pulse
        % (sometimes a few samples longer for larger flip pulses)
        Nset = length(rftemp{1});
        dur_set = Nset * SS_TS;
        rf = zeros(Nset,Nt); g = rf;
        for t = 1:Nt
            rf(1:Nset,t) = rftemp{t}(1:Nset);
            g(1:Nset,t) = gtemp{1}(1:Nset);
        end
        
        fprintf(1, 'Fs: %6.1f, Maximums - B1: %5.3fG Power: %5.3e G^2 ms Dur: %4.1fms\n', fs_best, b1_set, pow_set, dur_set*1e3);
    else
        fprintf(1, 'Fs: %6.1f *** No Soln ***\n', fs_best);
    end
else
    fprintf(1, 'Iterating on spectral sampling frequency to reduce B1\n');
    
    % Keep iterating on pulse design until B1 requirement is
    % met with highest sampling rate possible
    %
    dur_best = inf;
    b1_best = inf;
    nbands = length(fs_bands)/2;
    ss_max_duration_orig = SS_MAX_DURATION;  ss_min_order_orig = SS_MIN_ORDER;
    
    
    for band = nbands:-1:1,
        % Try each band
        %
        fs_bot = fs_bands(band*2-1);
        fs_top = fs_bands(band*2);
        d_idx = floor((fs_top - fs_bot + 1)/(SS_NUM_LOBE_ITERS-1));
        d_idx = max(1,d_idx);
        niter = ceil((fs_top-fs_bot+1)/d_idx);
        for iter = niter:-1:1,
            idx = fs_bot + (iter-1)*d_idx;
            if iter == niter,
                idx = fs_top;
            end;
            fs = fs_test(idx);
            
            % Call design script for SS that calculates rf, g that meets
            % sampling requirements
            %
            b1_set = 0;  dur_set = 0; pow_set = 0; pw_set = 0; fs_fail = 0;
            SS_MAX_DURATION = ss_max_duration_orig;  SS_MIN_ORDER =  ss_min_order_orig;
            
            for t = 1:Nt
                
                if strfind(ss_type, 'Flyback')
                    [rftemp{t}, gtemp{t}] = ss_flyback(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                else
                    [rftemp{t}, gtemp{t}] = ss_ep(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                end;
                
                if ~isempty(rftemp{t}),
                    pw_set = max([length(rftemp{t}), pw_set]);
                else
                    fs_fail = 1;
                end;
            end
            
            if fs_fail
                fprintf(1, 'Band: %d/%d Iter: %d/%d Fs: %6.1f *** No Soln ***\n', ...
                    nbands-band+1, nbands, niter-iter+1, niter, fs);
                continue;
            end;
            
            % fix pulse duration for all RF pulses
            ss_max_duration_orig = SS_MAX_DURATION;  ss_min_order_orig = SS_MIN_ORDER;
            SS_MAX_DURATION = pw_set*SS_TS;  SS_MIN_ORDER = 0;
            for t = 1:Nt
                
                if strfind(ss_type, 'Flyback')
                    [rftemp{t}, gtemp{t}] = ss_flyback(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                else
                    [rftemp{t}, gtemp{t}] = ss_ep(ang(t), z_thk, z_tb, z_d(:,t), f, a(:,t), d(:,t), fs, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                end;
                
                if ~isempty(rftemp{t}),
                    b1_set = max([abs(rftemp{t}),b1_set]);
                    dur_set = max([length(rftemp{t}) * SS_TS, dur_set]);
                    pow_set = max([sum(abs(rftemp{t}).^2) * SS_TS, pow_set]);
                else
                    fs_fail = 1;
                end;
                
            end
            
            if fs_fail
                fprintf(1, 'Band: %d/%d Iter: %d/%d Fs: %6.1f *** No Soln ***\n', ...
                    nbands-band+1, nbands, niter-iter+1, niter, fs);
                continue;
            end;
            
            nsolutions = nsolutions+1;
            
            
            % fix rephasing gradient area to be same as initial pulse
            % (sometimes a few samples longer for larger flip pulses)
            Nset = length(rftemp{1});
            dur_set = Nset * SS_TS;
            rf = zeros( Nset,Nt); g = rf;
            for t = 1:Nt
                rf( 1:Nset,t) = rftemp{t}(1:Nset);
                g(1:Nset,t) = gtemp{1}(1:Nset);
            end
            
            
            rfall{nsolutions} = rf;
            gall{nsolutions} = g;
            fsall(nsolutions) = fs;
            
            infoall{nsolutions} = sprintf('Fs: %6.1f maximums- B1: %5.3fG Power: %5.3e G^2 ms Dur: %4.1fms', ...
                fs, b1_set, pow_set, dur_set*1e3);
            fprintf(1, 'Band: %d/%d Iter: %d/%d %s\n', nbands-band+1, nbands, niter-iter+1, niter, infoall{nsolutions});
            
            
            if ((b1_set <= SS_MAX_B1) && (dur_set < dur_best))
                b1_best = b1_set;
                dur_best = dur_set;
                
                Ibest = nsolutions;
            elseif ((b1_best > SS_MAX_B1) && (b1_set < b1_best))
                b1_best = b1_set;
                dur_best = dur_set;
                
                Ibest = nsolutions;
            end
        end;
    end;
    
    % Choose desired solution
    %
    if nsolutions > 0
        fprintf(1, '\n');
        fprintf(1,'Solution(s) exists!\n');
        for n=1:nsolutions
            fprintf(1, '%d) %s\n', n, infoall{n});
        end
        if nsolutions > 1
            Isolution = input('Which pulses would you like to use? (leave empty for shortest pulse) ');
            if isempty(Isolution) || Isolution < 1 || Isolution > nsolutions
                Isolution = Ibest;
            end
            fprintf(1, 'Returning %s\n', infoall{Isolution});
        else
            Isolution = 1;
        end
        rf = rfall{Isolution};
        g = gall{Isolution};
        fs_best = fsall(Isolution);
    end
    
end;



if isempty(rf)
    error('No solution found! Try reducing bandwidths, increasing ripple, increasing max duration, increasing slice thickness...')
end;

if (no_plot==0)
    % Test RF
    %
    fmid = (f(1:2:end) + f(2:2:end))/2;
    bw = [min([f,-fs_best/2]) max([f,fs_best/2])];
    for t = [1, round(Nt/2) Nt]
        [f_plot,z_plot,m_plot] = ss_plot(g(:,t),rf(:,t),SS_TS, ptype,z_thk*2,bw,...
            SS_GAMMA, fmid);
        set(gcf, 'Name', ['Pulse #' int2str(t)])
    end
end


