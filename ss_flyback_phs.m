function [rf,g,isodelay] = ss_flyback_phs(ang, z_thk, z_tb, z_d, f, a, d, fs, ptype, z_ftype, ...
			     s_ftype, ss_type, f_off, dbg)
% SS_FLYBACK - Calculate flyback-type SS pulse
%   

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
% $Header: /home/adam/cvsroot/src/ss/ss_flyback_phs.m,v 1.4 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ss_globals;
    
    if strfind(ss_type,'Half')
	sym_flag = 1;
    else
	sym_flag = 0;
    end;

    % Check specification
    %
    [f_a, a_a, d_a, f_off] = ss_alias(f,a,d,f_off,fs,sym_flag);
    if (isempty(f_a))
	error('Strange: this frequency should be ok');
    end;
    if (dbg >= 2)
	ss_band_plot_phs(f_a, a_a, d_a, f_off, fs, min(f), max(f),sym_flag);
    end;
    
    % Calculate cycles/cm required
    %
    kz_max = z_tb / z_thk;		% cycles/cm
    kz_area = kz_max / SS_GAMMA;	% G/cm * s

    nsamp = round(1/(fs*SS_TS));
    [gpos, gneg, g1, g2, g3] = grad_ss(kz_area, nsamp, SS_VERSE_FRAC, SS_MXG, ...
				       SS_MXS, SS_TS, SS_EQUAL_LOBES);
    ng1 = length(g1);
    ng2 = length(g2);
    ng3 = length(g3);
    
    % Determine max order that can be supported 
    %
    t_poslobe = length(gpos) * SS_TS;
    t_lobe = length([gpos gneg]) * SS_TS;
    max_lobe = floor((SS_MAX_DURATION - t_poslobe) / t_lobe) + 1;
    
    % Prepare amplitude description that is consistent 
    % with other fir design calls
    a_dup = zeros(size(f_a));
    a_dup(1:2:end) = a_a;
    a_dup(2:2:end) = a_a;

    % Call fir filter design based on spectral factorization 
    % and convex optimization
    %
    if SS_MIN_ORDER, 
	use_max = 0;
    else
	use_max = 1;
    end;
    
    switch (s_ftype)
      case 'min'
        %      [s_b,status] = fir_minphase_power(max_lobe, f_a, a_dup, d_a,
        %      use_max, dbg);
        error('min phase spectral filter not supported for complex band specification');
      case 'max'
        %      [s_b,status] = fir_minphase_power(max_lobe, f_a, a_dup, d_a, use_max, ...
        %                                        dbg);
        %      s_b = conj(s_b(end:-1:1));
        error('max phase spectral filter not supported for complex band specification');
      case 'lin'
        if use_max, 
            % Make sure to use n_odd  - WHY?
            %
            %if bitget(max_lobe,1) == 0, 
            % max_lobe = max_lobe+1;
            % end;
            [s_b,status] = fir_qprog_phs(max_lobe, f_a, a_dup, d_a, [], dbg);
        else
            odd_or_even = 0;
            [s_b,status] = fir_min_order_qprog_phs(max_lobe, f_a, a_dup, d_a, odd_or_even, dbg);
        end
      case 'min_power'
        if use_max, 
            % Make sure to use n_odd -- WHY?
            %
            %            if bitget(max_lobe,1) == 0, 
            %                max_lobe = max_lobe+1;
            %            end;
            
            % Iterate over attempting different phase profiles for spectral
            % filter, choosing the one that minimizes the peak power
            %
            [s_b, a_dup_best, status] = fir_min_power_phs(max_lobe, f_a, a_dup, d_a, [], dbg);
        else
            error('Don''t handle case of minimum order, minimum peak power');
        end
        a_a = a_dup_best(1:2:end);
    end
    if strcmp(status, 'Solved');
	nlobe = length(s_b);
	if dbg >= 2,
	    ss_band_plot_phs(f_a, a_a, d_a, f_off, fs, min(f), max(f),sym_flag);

	    figure;
	    cplot(s_b);
	    title(sprintf('Filter Taps - Fs: %6.1f', fs));
	    xlabel('Tap');
	    ylabel('Amplitude');
	    drawnow;
	    pause(1);
	end;
	if dbg >= 1, 
	    figure;
	    nf = 512;
	    h = fftf(s_b, nf);
	    freqs = fs * [-nf/2:nf/2-1]/nf;
	    plot(freqs,abs(h));
	    title(sprintf('Spectral response - Fs: %6.1f', fs));
	    xlabel('Frequency');
	    ylabel('Magnitude');
	    grid;
	    drawnow;
	    pause(1);
	end;
	
	% Get Z RF pulse
	%
	if (dbg)
	    fprintf(1,'Getting Z RF pulse\n');
	end;
	z_np = length(g2);
	z_b = dzbeta(z_np, z_tb, 'st', z_ftype, ...
		      z_d(1), z_d(2));

	% Correct for non-linear effects with SLR if desired
	%
	if (SS_SLR_FLAG == 1), 
	    % Calculate excitation profile assuming in small-tip regime
	    %
	    oversamp = 4;
	    nZ = oversamp * z_np;
	    nZ2 = 2^ceil(log2(nZ));
	    Z_b = fftf(z_b, nZ2);  % column transform, unit magnitude
	    
	    if SS_SPECT_CORRECT_FLAG, 
		% Interpolate spectral filter on a grid that's equal to 
		% the number of time samples of the gradient --- do this 
		% partly before calling b2rf and partly afterwards
		%
		Nper = length([gpos gneg]);
		oversamp_slr = 16;
		Ntotal  = Nper * (length(s_b));
%		s_bi = 1/oversamp_slr * interpft([s_b], oversamp_slr * ...
%						 (length(s_b)));
		off = floor(z_np/2)/Nper;
%		Ntotal  = Nper * (length(s_b));
%		s_bi = 1/oversamp_slr * interpft([s_b], oversamp_slr * (length(s_b)));
		%		s_bi = length(s_b)/Ntotal * interpft(s_b,Ntotal); % Scale to keep
		% tranform consistent
	    
		% For each Z position: 
		%     - calculate nominal scaling of s_b
		%     - determine spectral taps through SLR
		%    
		bsf = sin(ang/2) * Z_b;
		s_rfm = zeros(length(bsf),length(s_b)*oversamp_slr);
		if 1
		    s_bi = spec_interp(s_b, oversamp_slr,-off,f_a, dbg);
		    for idx = 1:nZ2, 
			for bidx = 1:oversamp_slr,
			    s_rfm(idx,bidx:oversamp_slr:end) = ...
				b2rf(bsf(idx) * s_bi(bidx:oversamp_slr:end));
			end;
		    end;
		else
		    s_bi = spec_interp2(s_b, oversamp_slr,-off);
		    for idx = 1:nZ2, 
			s_rfm(idx,:) = oversamp_slr * b2rf(bsf(idx) * s_bi);
		    end;
		end;

		% Now calculate new z beta polynomial by inverse FFT, then use
		% SLR to get RF
		%
		z_bm = fftr(sin(s_rfm/2), z_np, 1); % Each row now scaled by
								   % tip-angle
		% Now do SLR in Z direction
		%
		if (dbg)
		    fprintf(1,'Doing SLR in Z...  \n');
		end;
		if 1 
		    z_rf = zeros(size(z_bm));
		    for idx = 1:size(z_bm,2), 
			z_rf(:,idx) = conj(b2rf(z_bm(:,idx)));
		    end;
		    
		    % Now raster scan for actual sampling positions
		    %
		    z_rfm = [];
		    for idx = 1:length(s_b)
			for zidx = 1:z_np, 
			    % idx_intrp = (idx*Nper - round(z_np/2) + zidx-1)/Ntotal*length(s_bi);
			    idx_intrp = 1 + ((idx-1)*Nper + zidx-1)/Ntotal*length(s_bi);
			    tmp_z_rf(zidx) = interp1([1:size(z_rf,2)], z_rf(zidx,:), ...
						     idx_intrp, 'spline');
			end;
			z_rfm = [z_rfm; tmp_z_rf(:).'];
		    end;
		else
		    % Choose time samples corresponding to trajectory
		    %
		    z_rfm = [];
		    for idx = 1:length(s_b), 
			for zidx = 1:z_np, 
			    %			idx_intrp = (idx*Nper - round(z_np/2) + zidx-1)/Ntotal*length(s_bi);
			    idx_intrp = 1 + ((idx-1)*Nper + zidx-1)/Ntotal*length(s_bi);
			    tmp_z_bm(zidx) = ...
				interp1([1:size(z_bm,2)], z_bm(zidx,:), ...
					idx_intrp, 'spline');
			end;
			tmp_z_rf = conj(b2rf(tmp_z_bm));
			z_rfm = [z_rfm; tmp_z_rf(:).'];
		    end;
		end;
        else % no spectral correction
		if (dbg)
		    fprintf(1,'Doing SLR in F...\n');
		end;
		s_rfm = [];
		bsf = sin(ang/2) * Z_b;
	    for idx = 1:nZ2, 
		    tmp_s_rf = b2rf(bsf(idx) * s_b); 
		    s_rfm = [s_rfm tmp_s_rf(:)];
		end;

		% Now calculate new z beta polynomial by inverse FFT, then use
		% SLR to get RF
		%
		z_bm = fftr(sin(s_rfm/2), z_np, 2); % Each row now scaled by
						    % tip-angle
		z_rfm = [];
		if (dbg)
		    fprintf(1,'Doing SLR in Z...  \n');
		end;
		for idx = 1:size(z_bm,1), 
		    tmp_z_rf = conj(b2rf(z_bm(idx,:)));
		    z_rfm = [z_rfm; tmp_z_rf(:).'];
		end;
	    end;

	    % Modulate rf to passband frequency BEFORE versing, then 
	    % modulate back. This will make sure that the slice profile
	    % shows no blurring at the passband.  In the case that 
	    % multiple passbands are defined, then the midpoint of  
	    % the first passband is used
	    %
	    pass_idx = find(a > 0, 1, 'first');
	    fpass = [f(pass_idx*2-1) f(pass_idx*2)];
	    fpass_mid = mean(fpass) - f_off;
	    
	    nlobe = length(s_b);
	    rfmod = exp(i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);

	    if SS_VERSE_B1
		if (dbg)
		    fprintf(1,'Versing RF with B1 minimization...  \n');
		end;

		% Additional VERSE with B1 restriction, maintaining duration
                if 1
                  % Verse based on largest RF amplitudes at each kz, then use same
                  % gradient lobe for each 
                  %
                  z_rfmax = max(abs(z_rfm));
                else 
                  % Verse largest RF pulse, then use same gradient lobe for
                  % each 
                  [b1max_sc,b1max_idx] = max(max(abs(z_rfm),[],2));
                  z_rfmax = z_rfm(b1max_idx,:);
                end
            
                z_rfvmax1 = ss_verse(g2, z_rfmax);

                %		z_rfmod = z_rfm(b1max_idx,:) .* rfmod;
		%z_rfvmod1 = ss_verse(g2, z_rfmod);
		[z_rfvmax, g2v] = ss_b1verse(g2, z_rfvmax1(:).', SS_MAX_B1, SS_MXG, ...
					    SS_MXS, SS_TS, SS_GAMMA, SS_SLEW_PENALTY, dbg);
		if (isempty(z_rfvmax))
		    % B1 condition cannot be met
		    rf = [];
		    g = [];
		    return;
        end

        if 0
            % experimental filtering:
            % cutoff at x*125 kHz
            b = firls(10, [0 .1 .2 1], [1 1 0 0]);
            g2vs = sum(g2v);
            g2v = filtfilt(b, 1, g2v); 
            g2v = g2v * g2vs / sum(g2v);
        end
        
		for idx = 1:nlobe, 
		    z_rfmod = z_rfm(idx,:) .* rfmod;
		    z_rfvmod = ss_verse(g2v, z_rfmod);
		    z_rfv = z_rfvmod(:).' .* conj(rfmod);
		    z_rfmv(idx,:) = z_rfv;
		end;

		% update gradient
		gpos = [g1, g2v, g3];
		if SS_EQUAL_LOBES
		    gneg = -gpos;
		end	
	    else
		if (dbg)
		    fprintf(1,'Versing RF...  \n');
		end;
		for idx = 1:nlobe, 
		    z_rfmod = z_rfm(idx,:) .* rfmod;
		    if SS_VERSE_FRAC == 0, 
			z_rfvmod = z_rfmod;
		    else
			z_rfvmod = ss_verse(g2, z_rfmod);
		    end;
		    z_rfv = z_rfvmod(:).' .* conj(rfmod);
		    z_rfmv(idx,:) = z_rfv;
		end;
	    end;
    else % no SLR
	    if SS_SPECT_CORRECT_FLAG, 
		% Spectral correction needs to be applied to unaliased bands, 
		% therefore the raw frequency spec needs to be passed through.
		%
%		s_b = [0; s_b; 0];  % play with adding extra taps to make 
		% spectral correction easier
		Nper = length([gpos gneg]);
		st_off = -floor(ng2/2);
		Noff = [st_off:st_off+ng2-1];	       
		bsf = sin(ang/2) * ones(size(Noff));	% This needs to be
                                                        % updated 
		if (dbg)
		    fprintf(1,'No SLR.. spectral correction...\n');
		end;
%		s_rfm = ss_spect_correct(s_b, bsf, Nper, Noff, (f-f_off)/(fs/2), ...
%					 ptype, 'Flyback', SS_SLR_FLAG, SS_SPECT_CORRECT_REGULARIZATION, dbg);
		s_rfm = ss_spect_correct(s_b, bsf, Nper, Noff, f, ...
					 ptype, 'Flyback', SS_SLR_FLAG, SS_SPECT_CORRECT_REGULARIZATION, dbg);
	    else
		s_rfm = ang * conj(s_b(:)) * ones(1,ng2); % Intentional conjugation
							  % here -- needed because of 
							  % possible asymmetric frequency
							  % response 
	    end

	    % Modulate rf to passband frequency BEFORE versing, then 
	    % modulate back. This will make sure that the slice profile
	    % shows no blurring at the passband.  In the case that 
	    % multiple passbands are defined, then the midpoint of  
	    % the first passband is used
	    %
	    pass_idx = find(a > 0, 1, 'first');
	    fpass = [f(pass_idx*2-1) f(pass_idx*2)];
	    fpass_mid = mean(fpass) - f_off;
	    z_bmod = z_b(:).' .* exp(i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);
	    
	    if SS_VERSE_B1
		if (dbg)
		    fprintf(1,'Versing RF with B1 minimization...  \n');
		end;
		% Additional VERSE with B1 restriction, maintaining duration
		% account for scaling by the largest possible spectral weightings at
		% each spatial sample
		b1max_sc = max(abs(s_rfm));
		
		z_bvmod1 = ss_verse(g2, z_bmod);
		[z_bvmod, g2v] = ss_b1verse(g2, z_bvmod1(:).', SS_MAX_B1 ./ b1max_sc, SS_MXG, ...
					  SS_MXS, SS_TS, SS_GAMMA, SS_SLEW_PENALTY, dbg);
		
		if (isempty(z_bvmod))
		    % B1 condition cannot be met
		    rf = [];
		    g = [];
		    return;
        end

        if 0
            % experimental filtering:
            % cutoff at x*125 kHz
            b = firls(10, [0 .1 .2 1], [1 1 0 0]);
            g2vs = sum(g2v);
            g2v = filtfilt(b, 1, g2v); 
            g2v = g2v * g2vs / sum(g2v);

            z_bvmod = ss_verse(g2v, z_bmod);
        end
        
		% update gradient
		gpos = [g1, g2v, g3];
		if SS_EQUAL_LOBES
		    gneg = -gpos;
		end	
	    else
		if (dbg)
		    fprintf(1,'Versing RF...  \n');
		end;
		if SS_VERSE_FRAC == 0, 
		    z_bvmod = z_bmod;
		else
		    z_bvmod = ss_verse(g2, z_bmod);
		end;
	    end;
	    
	    % Modulate back
	    %
	    z_bv = z_bvmod(:).' .* exp(-i*2*pi*[0:(z_np-1)]*SS_TS* ...
				       fpass_mid);
	    
	    % Build RF matrix
	    %
	    nlobe = length(s_b);
	    z_rfmv = s_rfm .* (ones(nlobe,1) * z_bv);
	end;
	
        % Calculate isodelay --- since we only have linear phase pulses for 
        % this routine, it's easy...
        %
        isodelay = ((nlobe/2) * length(gpos) + (nlobe-1)/2*length(gneg)) ...
            * SS_TS;
            
	% Compile g, RF
	%
	rf = []; g = [];
	nneg = length(gneg);
	for idx = 1:nlobe,
	    rf_lobe = z_rfmv(idx,:);
	    rf = [rf zeros(1,ng1) rf_lobe zeros(1,ng3)];
	    if idx < nlobe,
		rf = [rf zeros(1,nneg)];
		g = [g gpos gneg];
	    else
		g = [g gpos];
	    end;
	end;

	% Offset RF 
	%
	nrf = length(rf);
	rf = rf .* exp(-i*2*pi*[0:nrf-1]*SS_TS*f_off);
	
	% Convert amplitude to Gauss
	% 
	rf = rf / (2 * pi * SS_GAMMA * SS_TS);
	
	% Calculate refocussing lobe
	%
	switch (ptype)
	 case{'ex', 'se'}
	  ;
	 otherwise
	  return;
	end;
	
	% Step 1 - get passbands
	%
	fmid = (f(1:2:end) + f(2:2:end))/2;
	idx_pass = find(a > 0);
	fpass = fmid(idx_pass);
	npass = length(fpass);
	
	% Step 2 - get spatial sample points
	%
	nz = 101;
	dz = z_thk / (nz-1);
	z = [-z_thk/2:dz:z_thk/2];
	
	% Step 3 - get spatial profile
	%
	gzrot = 2 * pi * SS_GAMMA * SS_TS * g;
	gfrot = 2 * pi * SS_TS * ones(size(g));
	rrot = 2 * pi * SS_GAMMA * SS_TS * rf;
	switch (ptype)
	 case {'ex'}
	  mxy = ab2ex(abr(rrot, gzrot + i*gfrot, z, fpass));
	 case {'se'}
	  mxy = ab2se(abr(rrot, gzrot + i*gfrot, z, fpass));
	 otherwise
	end;
	
	% Step 4 - find best fit to phase ramp
	%
	zpass = z(:) * ones(1,npass);
	zpass = zpass(:);
	mxy_phs = unwrap(angle(mxy));
	mxy_phs_mid = mxy_phs((nz+1)/2,:);
	mxy_phs = mxy_phs - (ones(nz,1) * mxy_phs_mid);
	mxy_phs = mxy_phs(:);
	
	p = polyfit(zpass,mxy_phs,1);
	slope = p(1);
	
	% Step 5 - get area of gradient required
	%
	cyc_per_cm = slope/(2*pi);
	m0 = cyc_per_cm / SS_GAMMA;
	gref = grad_mintrap(m0, SS_MXG, SS_MXS, SS_TS);
	
	rf = [rf zeros(size(gref))];
	g = [g gref];
        isodelay = isodelay + length(gref) * SS_TS;
    else
	% No solution for this frequency
	%
	rf = [];
	g = [];
        isodelay = [];
    end;
    
