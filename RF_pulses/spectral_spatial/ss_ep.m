function [rf,g] = ss_ep(ang, z_thk, z_tb, z_d, f, a, d, fs, ptype, z_ftype, ...
			     s_ftype, ss_type, f_off, dbg)
% SS_EP - Calculate echo-planar-type SS pulse
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
% $Header: /home/adam/cvsroot/src/ss/ss_ep.m,v 1.16 2013/08/15 16:02:33 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ss_globals;
    
    if strfind(ss_type,'Half')
	sym_flag = 1;
    else
	sym_flag = 0;
    end;
          
    [f_a, a_a, d_a, f_off] = ss_alias(f,a,d,f_off,fs,sym_flag);
    if (isempty(f_a))
	error('Strange: this frequency should be ok');
    end;
    if (dbg >= 3)
	ss_band_plot(f_a, a_a, d_a, f_off, fs, min(f), max(f),sym_flag);
    end;
    
    % Calculate cycles/cm required
    %
    kz_max = z_tb / z_thk;		% cycles/cm
    kz_area = kz_max / SS_GAMMA;	% G/cm * s

    % Calculate SS bipolars
    %
    nsamp = round(2/(fs*SS_TS));
    [gpos, gneg, g1, g2, g3] = grad_ss(kz_area, nsamp, SS_VERSE_FRAC, SS_MXG, ...
				       SS_MXS, SS_TS, 1);
    ng1 = length(g1);
    ng2 = length(g2);
    ng3 = length(g3);
    
    % Determine max order that can be supported 
    %
    t_lobe = length(gpos) * SS_TS;
    max_lobe = floor(SS_MAX_DURATION / t_lobe);

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
      [s_b,status] = fir_minphase_power(max_lobe, f_a, a_dup, d_a, use_max, dbg);
     case 'max'
      [s_b,status] = fir_minphase_power(max_lobe, f_a, a_dup, d_a, use_max, ...
                                        dbg);
      s_b = conj(s_b(end:-1:1));
     case 'lin'
      if use_max, 
	  % Make sure to use n_odd 
	  %
	  if bitget(max_lobe,1) == 0, 
	      max_lobe = max_lobe+1;
	  end;
	  [s_b,status] = fir_qprog(max_lobe, f_a, a_dup, d_a, dbg);
      else
	  odd_or_even = 0;
	  [s_b,status] = fir_min_order_qprog(max_lobe, f_a, a_dup, d_a, odd_or_even, dbg);
      end;
    end;
    if strcmp(status, 'Solved');
	nlobe = length(s_b);
	if dbg >= 2,
	    ss_band_plot(f_a, a_a, d_a, f_off, fs, min(f), max(f),sym_flag);

	    figure;
	    cplot(s_b);
	    title(sprintf('Filter Taps - Fs: %6.1f', fs));
	    xlabel('Tap');
	    ylabel('Amplitude');
	    drawnow;
	    pause(1);
	    
	    figure;
	    nf = 512;
	    h = fftshift(fft(s_b, nf));
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
            fprintf(1,'Spectral Correction not supported for EP pulses with SLR\n');
        end
        if 0 %SS_SPECT_CORRECT_FLAG
        % Currently not working
    	% Interpolate spectral filter on a grid that's equal to 
		% the number of time samples of the gradient --- do this 
		% partly before calling b2rf and partly afterwards
		%
		Nper = length(gpos);
		oversamp_slr = 16;  % Increasing this helps...
		Ntotal  = Nper * (length(s_b));
		off = (z_np/2)/Nper;
	    
		% For each Z position: 
		%     - calculate nominal scaling of s_b
		%     - determine spectral taps through SLR
		%    
		bsf = sin(ang/2) * Z_b;

        
        % non-uniform spectral interpolation
        % The sampling of the output is different than spec_interp
        [s_bi_nu, t_nu] = spec_interp_nonuniform(s_b, oversamp_slr, off,f_a, dbg);
            
        s_rfm_nu = zeros(length(bsf),length(s_bi_nu));
        % Non-uniform SLR would be good here...
        for zidx = 1:nZ2
            if 0
                % small-tip in spectral direction
                s_rfm_nu(zidx,:) = ang*Z_b(zidx) * conj(s_bi_nu); % conj?
                %s_rfm_nu(zidx,:) = 2*asin(abs(bsf(zidx))) * conj(s_bi_nu); % conj?
                % s_rfm_nu(zidx,:) = ang*abs(Z_b(zidx)) * conj(s_bi_nu);
            else
                % SLR on uniformly spaced samples
                [t_nu_sort, I_nu_sort] = sort(t_nu);
                t_nu = t_nu_sort;
                s_bi_nu = s_bi_nu(I_nu_sort);
                for idx = 1:oversamp_slr
                    s_rfm_nu(zidx,idx:oversamp_slr:end) = ...
                        b2rf(bsf(zidx) * s_bi_nu(idx:oversamp_slr:end));
                end
                
            end
        end

        % Now calculate new z beta polynomial by inverse FFT
		z_bm_nu = fftr(sin(s_rfm_nu/2), z_np, 1); % Each row now scaled by
								   % tip-angle

        % SLR in Z to get RF
	    z_rf_nu = zeros(size(z_bm_nu));
	    for idx = 1:size(z_bm_nu,2), 
			z_rf_nu(:,idx) = conj(b2rf(z_bm_nu(:,idx)));
	    end;
		    
	    % Now raster scan for actual sampling positions
	    z_rfm = [];
        
        for idx = 1:length(s_b)
			for zidx = 1:z_np,
                zi_off = (2*zidx - z_np-1) / (z_np-1) * off;
                if (rem(idx,2) == 1)
                    idx_intrp = (idx-1) + zi_off;
                else
                    idx_intrp = (idx-1) - zi_off;
                end
                tmp_z_rf(zidx) = interp1(t_nu, ...
                    z_rf_nu(zidx,:), idx_intrp, 'spline');
        	end;
            z_rfm = [z_rfm; tmp_z_rf(:).'];
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
	    end; %SS_SPECT_CORRECT_FLAG

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
            %            
            if 1
              % Verse based on largest RF amplitudes at each kz, then use same
              % gradient lobe for each 
              %
              z_rfmax = max(abs(z_rfm));
            else 
              % verse based on largest RF subpulse
              [b1max_sc,b1max_idx] = max(max(abs(z_rfm),[],2));
              z_rfmax = z_rfm(b1max_idx,:);
            end
            
            z_rfvmax1 = ss_verse(g2, z_rfmax);
            [z_rfvmax, g2v] = ss_b1verse(g2, z_rfvmax1(:).', SS_MAX_B1, SS_MXG, ...
                            SS_MXS, SS_TS, SS_GAMMA, SS_SLEW_PENALTY, dbg);
            

            if (isempty(z_rfvmax))
                % B1 condition cannot be met
                rf = [];
                g = [];
                return;
            end

            for idx = 1:nlobe, 
              if (rem(idx,2) == 1)
                z_rfmod = z_rfm(idx,:) .* rfmod;
                z_rfvmod = ss_verse(g2v, z_rfmod);
              else
                z_rfmod = z_rfm(idx,end:-1:1) .* rfmod;
                z_rfvmod = ss_verse(g2v(end:-1:1), z_rfmod);
              end
              z_rfv = z_rfvmod(:).' .* conj(rfmod);
              z_rfmv(idx,:) = z_rfv;
            end;

            % update gradient
            gpos = [g1, g2v, g3];
            gneg = [-g1, -g2v(end:-1:1), -g3];
	    else
            if (dbg)
                fprintf(1,'Versing RF...  \n');
            end;
            for idx = 1:nlobe, 
                if (rem(idx,2) == 1)
                  z_rfmod = z_rfm(idx,:) .* rfmod;
                else
                  z_rfmod = z_rfm(idx,end:-1:1) .* rfmod;
                end
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
        % Correct spectral filter for non-uniform sampling 
        % in spectral dimension
        %
        Nper = length([gpos]);
        Noff = -(ng2-1)/2:((ng2-1)/2);	       
        bsf = sin(ang/2) * ones(size(Noff));
        if SS_SPECT_CORRECT_FLAG, 
        %	    s_b = [0; s_b];	% Can pad to keep from blowing up

            s_rfm = ss_spect_correct(s_b, bsf, Nper, Noff, (f-f_off)/(fs/2), ...
                         ptype, 'EP', SS_SLR_FLAG, dbg);
            s_rfm = conj(s_rfm);
        else
            s_rfm = ang * conj(s_b(:)) * ones(1,ng2); % Intentional conjugation
                                  % here -- needed because of 
                                  % possible asymmetric frequency
                                  % response 
        end;


        % Modulate rf to passband frequency BEFORE versing, then 
        % modulate back. This will make sure that the slice profile
        % shows no blurring at the passband.  In the case that 
        % multiple passbands are defined, then the midpoint of  
        % the first passband is used.  
        %
        % NOTE:  With an EP trajectory, we have to play the RF 
        % out in reverse.  As as result, we should make sure to 
        % verse the RF pulse in both directions separately in case
        % it is not symmetric. 
        %
        % VERSE'ing in both forward and reverse directions has the effect of
        % introducing significant distortions into the other bands, which is
        % particularly bad in the stop bands. Not sure if there is a way around
        % this, but only VERSEing in one direction and then using a reversed
        % pulse seems to reduce this problem.

        pass_idx = find(a > 0, 1, 'first');
        fpass = [f(pass_idx*2-1) f(pass_idx*2)];
        fpass_mid = mean(fpass) - f_off;
        fpass_mid = 0;
        z_bmod_for = z_b(:).' .* exp(i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);
        z_b_rev = z_b(end:-1:1);
        z_bmod_rev = z_b_rev(:).' .* exp(i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);


        if SS_VERSE_B1
            if (dbg)
                fprintf(1,'Versing RF with B1 minimization...  \n');
            end;
            % Additional VERSE with B1 restriction, maintaining duration
            % account for scaling by the largest possible spectral weightings at
            % each spatial sample
            b1max_sc = max(abs(s_rfm), [], 1);

            z_bvmod1 = ss_verse(g2, z_bmod_for);
            [z_bvmod_for, g2v_for] = ss_b1verse(g2, z_bvmod1(:).', SS_MAX_B1 ./ b1max_sc, SS_MXG, ...
                          SS_MXS, SS_TS, SS_GAMMA, SS_SLEW_PENALTY, dbg);

            % not sure this is necessary, adds nice symmetry though
            z_bvmod1 = ss_verse(g2, z_bmod_rev);
            [z_bvmod_rev, g2v_rev] = ss_b1verse(g2, z_bvmod1(:).', SS_MAX_B1 ./ b1max_sc(end:-1:1), SS_MXG, ...
                          SS_MXS, SS_TS, SS_GAMMA, SS_SLEW_PENALTY, dbg);

            if (isempty(z_bvmod_for) || isempty(z_bvmod_rev))
                % B1 condition cannot be met
                rf = [];
                g = [];
                return;
            end

            g2v = min(g2v_for, g2v_rev(end:-1:1));

            z_bvmod_for = ss_verse(g2v, z_bmod_for);
            z_bvmod_rev = ss_verse(g2v(end:-1:1), z_bmod_rev);

            z_bv_for = z_bvmod_for(:).' .* exp(-i*2*pi*[0:(z_np-1)]*SS_TS* ...
                           fpass_mid);
            z_bv_rev = z_bvmod_rev(:).' .* exp(-i*2*pi*[0:(z_np-1)]*SS_TS* ...
                           fpass_mid);
            z_bv_rev = z_bv_for(end:-1:1);


            % update gradient
            gpos = [g1, g2v, g3];
            gneg = [-g1, -g2v(end:-1:1), -g3];
        else


            % Do forward RF pulse 
            %
            z_bvmod_for = ss_verse(g2, z_bmod_for);	

            % Modulate back to 0
            %
            z_bv_for = z_bvmod_for(:).' .* exp(-i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);

            % Do reverse RF pulse
            %
            z_bvmod_rev = ss_verse(g2, z_bmod_rev);
            z_bv_rev = z_bvmod_rev(:).' .* exp(-i*2*pi*[0:(z_np-1)]*SS_TS*fpass_mid);

            % More accurate to just reverse RF pulse
            % reverse VERSE'ing doesn't work well with asymmetric RF pulses for
            % some reason
            z_bv_rev = z_bv_for(end:-1:1);
        end	    

	    % Build RF matrix
	    %
	    nlobe = length(s_b);
            z_rfmv = zeros(nlobe, z_np);
            for idx = 1:nlobe
              if rem(idx,2) == 1
		z_rfmv(idx,:) = s_rfm(idx,:) .* z_bv_for;
              else
		z_rfmv(idx,:) = s_rfm(idx,end:-1:1) .* z_bv_rev;
              end;
            end
            
    end;  % SLR
	
    
	% Compile g, RF
	%
	rf = []; g = [];
	for idx = 1:nlobe,
          rf_lobe = z_rfmv(idx,:);
          rf = [rf zeros(1,ng1) rf_lobe zeros(1,ng3)];
          if rem(idx,2) == 1
            g = [g gpos];
          else
            g = [g gneg];
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
    else
	% No solution for this frequency
	%
	rf = [];
	g = [];
    end;
    
