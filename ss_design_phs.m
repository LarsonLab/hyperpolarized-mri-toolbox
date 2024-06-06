function [g, rf, fs_best, z_plot, f_plot, m_plot, isodelay] = ...
    ss_design_phs(z_thk, z_tb, z_de, f, a_angs, de, a_phs, d_phs, ... 
	      ptype, z_ftype, s_ftype, ss_type, ...
	      f_off, dbg);
    
					
% SS_DESIGN - Design spectral-spatial pulse
%   
% [g, rf, fs_best, z_plot, f_plot, m_plot] = ...
%     ss_design(z_thk, z_tb, z_de, f, a_angs, de,...
%               ptype, z_ftype, s_ftype, ss_type, ...
%               f_off, dbg)
%
% INPUTS
%   z_thk - slice thickness (cm)
%   z_tb - spatial time-bandswidth
%   z_de - spatial ripples, [pass_ripple, stop_ripple]
%   f - spectral band edge specification (Hz)
%   a_ang - spectral band flip angle specification (radians)
%   de - spectral band ripples - inphase
%   a_phs - spectral band phase specification (radians)
%   d_phs - spectral band ripples - phase - (radians)
%   ptype - spatial pulse type: 'ex' (default), 'se', 'sat', 'inv'
%   z_ftype - spatial filter type: 'ms', 'ls', 'pm' (default), 'min', 'max'
%   s_ftype - spectral filter type: 'min' (default), 'max', 'lin'
%   ss_type - spectral-spatial type: 'Flyback Whole' (default),
%       'Flyback Half', 'EP Whole', 'EP Half',
%       'EP Whole Opp-Null', 'EP Half Opp-Null'
%   f_off - center frequency (empty to let ss_design choose)
%   dbg - print debug messages: 0-none (default), 1-some, 2-more
%
% OUTPUTS
%   g - gradient (G/cm)
%   rf - RF (G)
%   fs_best - spectral sampling frequency (Hz)
%   z_plot, f_plot, m_plot - ss_plot() outputs, ranges and values in figure
%   isodelay - time from effective tipdown to end of pulse
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
% $Header: /home/adam/cvsroot/src/ss/ss_design_phs.m,v 1.5 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check all inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    if (nargin < 8), 
	error(['Usage: ss_design(z_thk, z_tb, z_d, f, a_angs, d, a_phs, d_quad' ...
	       'ptype, z_fttype, s_ftype, ss_type, foff)']); 
    end;
    
    % Check ptype
    %
    if (nargin < 9) || isempty(ptype), 
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
    if (nargin < 10) || isempty(z_ftype), 
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
    if (nargin < 11) || isempty(s_ftype), 
	s_ftype = 'min'; 
    else
	switch s_ftype, 
	 case {'min', 'max', 'lin', 'min_power'}
	 otherwise 
	  error(sprintf(['Spectral filter type (s_ftype) of: %s not' ...
			 ' recognized'], s_ftype));
	end;
    end;
    
    % Check ss_type
    %
    if (nargin < 12) || isempty(ss_type), 
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
    if (nargin < 13) || isempty(f_off), 
	f_off = [];
    end;
    
    % Check dbg
    %
    if (nargin < 14) || isempty(dbg), 
	dbg = 0;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize globals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ss_globals;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert "a_angs" to Beta polynomial "a"
    % Convert effective ripples to polynomial ripples 
    % depending on pulse type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ang = max(a_angs);
    
    if SS_SLR_FLAG, 
	a = sin(a_angs/2)/sin(ang/2);
	%a = asin(sin(ang)*a)/ang;
    else
	a = a_angs/ang;
    end
    
    z_d = rf_ripple(z_de, [1 0], ang, ptype);
    [d,a,ang] = rf_ripple(de, a, ang, ptype);

    % Update d, a to be complex to include phase specification
    % 
    d = d .* exp(i*d_phs);
    a = a .* exp(i*a_phs);
    
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
    for idx = 1:SS_NUM_FS_TEST,
	fs_test(idx) = fs_min + (idx-1)*df;
	% Make sure fs_cur is integer number of samples
	%
	nsamp = ceil(1/(fs_test(idx)*SS_TS));
	fs_test(idx) = 1/(nsamp*SS_TS);

	% Test aliasing of frequencies into operating BW
	%
	[f_a, a_a, d_a, fo] = ss_alias(f,a,d,f_off,fs_test(idx)*fs_scale,sym_flag);
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
	    [f_a, a_a, d_a, fo] = ss_alias(f,a,d,f_off,fs_over*fs_scale,sym_flag);
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

	if strfind(ss_type, 'Flyback')
	    [rf, g, isodelay] = ss_flyback_phs(ang, z_thk, z_tb, z_d, f, a, d, fs_best, ...
					   ptype, z_ftype, s_ftype, ss_type, ...
					   f_off, dbg);
	else
	    [rf, g, isodelay] = ss_ep_phs(ang, z_thk, z_tb, z_d, f, a, d, fs_best, ...
				      ptype, z_ftype, s_ftype, ss_type, ...
				      f_off, dbg);
	end;
	    
	
	if ~isempty(rf),
	    b1_best = max(abs(rf));
	    dur_best = length(rf) * SS_TS;
        pow_best = sum(abs(rf).^2) * SS_TS;
        nsolutions = 1;
    	fprintf(1,'Solution(s) exists!\n');
    	fprintf(1, 'Fs: %6.1f B1: %5.3fG Power: %5.3e G^2 ms Dur: %4.1fms\n', fs_best, b1_best, pow_best, dur_best*1e3); 
    else
        fprintf(1, 'Fs: %6.1f *** No Soln ***\n', fs_best);

	end;
    else
	fprintf(1, 'Iterating on spectral sampling frequency to reduce B1\n');

	% Keep iterating on pulse design until B1 requirement is 
	% met with highest sampling rate possible
	%
	dur_best = inf;
	b1_best = inf;
	nbands = length(fs_bands)/2;
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
		if strfind(ss_type, 'Flyback')
		    [rf, g, isodelay] = ss_flyback_phs(ang, z_thk, z_tb, z_d, f, a, d, fs, ...
					 ptype, z_ftype, s_ftype, ss_type, ...
					 f_off, dbg);
		else
		    [rf, g, isodelay] = ss_ep_phs(ang, z_thk, z_tb, z_d, f, a, d, fs, ...
				    ptype, z_ftype, s_ftype, ss_type, ...
				    f_off, dbg);
		end;
		    
		if isempty(rf), 
		    fprintf(1, 'Band: %d/%d Iter: %d/%d Fs: %6.1f *** No Soln ***\n', ...
			    nbands-band+1, nbands, niter-iter+1, niter, fs);
		    continue;
		end;
		
        nsolutions = nsolutions+1;
        
        rfall{nsolutions} = rf;
        gall{nsolutions} = g;
        fsall(nsolutions) = fs;
        isodelayall(nsolutions) = isodelay;
        
		b1 = max(abs(rf));
		dur = length(rf) * SS_TS;
        pow = sum(abs(rf).^2) * SS_TS;
        
        infoall{nsolutions} = sprintf('Fs: %6.1f B1: %5.3fG Power: %5.3e G^2 ms Dur: %4.1fms', ...
			fs, b1, pow, dur*1e3);
		fprintf(1, 'Band: %d/%d Iter: %d/%d %s\n', nbands-band+1, nbands, niter-iter+1, niter, infoall{nsolutions}); 
        
        
         if ((b1 <= SS_MAX_B1) && (dur < dur_best))
 		    b1_best = b1;
 		    dur_best = dur;

            Ibest = nsolutions; 		
 		elseif ((b1_best > SS_MAX_B1) && (b1 < b1_best))
 		    b1_best = b1;
 		    dur_best = dur;

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
            Isolution = input('Which pulse would you like to use? (leave empty for shortest pulse) ');
            if isempty(Isolution) || Isolution < 1 || Isolution > nsolutions
                Isolution = Ibest;
            end
            fprintf(1, 'Returning %s\n', infoall{Isolution});
        else
            Isolution = 1;
        end
        rf = rfall{Isolution};
        g = gall{Isolution};
        isodelay = isodelayall(Isolution);
        fs_best = fsall(Isolution);
    end
    
    end;
      
    
	
    if isempty(rf)
        fprintf(1,'No solution found! Trying to increase band ripples to determine limiting\n frequency specifications...\n')
    
        orig_min_order = SS_MIN_ORDER;
        SS_MIN_ORDER = 0;
      	fs_top = fs_bands(end);
        fs_best = fs_test(fs_top);

        d_max = 1*ones(size(d));
        d_min = abs(d);
        tol_factor = 4;
        
        for Id = 1:length(d)
            % trying bisection search with increased ripple in each band to
            % determine which portion of frequency spec is limiting
            d_test = d;
            
            while (d_max(Id) - d_min(Id)) > abs(d(Id))/tol_factor
                d_test(Id) = (d_max(Id) + d_min(Id))/2 * exp(i*angle(d(Id)));
                
                if strfind(ss_type, 'Flyback')
                    [rf, g, isodelay] = ss_flyback_phs(ang, z_thk, z_tb, z_d, f, a, d_test, fs_best, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                else
                    [rf, g, isodelay] = ss_ep_phs(ang, z_thk, z_tb, z_d, f, a, d_test, fs_best, ...
                        ptype, z_ftype, s_ftype, ss_type, ...
                        f_off, dbg);
                end;
                
                % update min/max ripple values
                if ~isempty(rf),
                    d_max(Id) = abs(d_test(Id));
                else
                    d_min(Id) = abs(d_test(Id));
                end;
            end
        end
        
        % find ripple values that create solutions
        Ifix = find(d_max < 1);
        if ~isempty(Ifix)
            [tempmin Imin] = min(d_max - abs(d));
            d_fix = d; d_fix(Imin) = d_max(Imin)*exp(i*angle(d(Imin)));
            fprintf(1,'Solution found by increasing ripples in bands %s\n', int2str(Ifix));
            fprintf(1,'Returning pulse with increased ripple in band %d (%.1f to %.1f Hz):\n', Imin, f(2*Imin-1), f(2*Imin));
           
            if strfind(ss_type, 'Flyback')
                [rf, g, isodelay] = ss_flyback_phs(ang, z_thk, z_tb, z_d, f, a, d_fix, fs_best, ...
                    ptype, z_ftype, s_ftype, ss_type, ...
                    f_off, dbg);
            else
                [rf, g, isodelay] = ss_ep_phs(ang, z_thk, z_tb, z_d, f, a, d_fix, fs_best, ...
                    ptype, z_ftype, s_ftype, ss_type, ...
                    f_off, dbg);
            end;
            	    
            b1_best = max(abs(rf));
            dur_best = length(rf) * SS_TS;
            pow_best = sum(abs(rf).^2) * SS_TS;
         	fprintf(1, 'Fs: %6.1f B1: %5.3fG Power: %5.3e G^2 ms Dur: %4.1fms\n', fs_best, b1_best, pow_best, dur_best*1e3); 
           
            fprintf(1,'\nPulse specs should be modified by reducing bandwidths or increasing ripple in bands %s\n', int2str(Ifix));
            fprintf(1,'Increasing the max pulse duration or slice thickness may also help\n');

        
        else
            error('No solution found! Try reducing bandwidths, increasing ripple, increasing max duration, increasing slice thickness...')
        end
        
        
        SS_MIN_ORDER = orig_min_order;
    end;
    
    % Test RF
    %
    fmid = (f(1:2:end) + f(2:2:end))/2;
    bw = [-fs_best/2 fs_best/2];
    [f_plot,z_plot,m_plot] = ss_plot(g,rf,SS_TS, ptype,z_thk*2,bw,...
        SS_GAMMA, fmid, isodelay);
    


   