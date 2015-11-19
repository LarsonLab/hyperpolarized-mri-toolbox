function rfm = ss_spect_correct(b, bsf, Nper, Noff, f, ptype, ss_type, slr, ...
				reg_factor, dbg)
% SS_SPECT_CORRECT - Correct spectral filter for irregular sampling
%   
% function rf = ss_spect_correct(b, bsf, Nper, Noff, f, ptype, ss_type, slr, reg_factor, dbg)
%    
% Inputs: 
%   b - spectral filter design, normalized so that passband has value of 1
%   bsf - beta polynomial scale factors (normally sin(ang/2))
%   Nper - period in samples between taps of b
%   Noff - vector of sample offsets from reference point (may be fractional)
%   f - Normalized frequency bands to try to correct - may be outside Nyquist (-1..1)
%   ptype - type of pulse, 'ex', 'se', 'sat, 'inv'
%   ss_type - 'Flyback' or 'EP'
%   slr - SLR flag
%   reg_factor - regularization factor to reduce peak amplitude increases
%      from matrix inversion in EP designs
%
% Outputs: 
%   rf - rf taps to give desired tip    
%   
% It is assumed that the first gradient lobe is a positive one, 
%
% If a "Flyback" ss_type is specified, then filter taps will all 
% be moving forward together by "n" units.  With minimum-phase filters
% it is best to use the spectral filter reference for the first tap
% then move forward.
%
% If a "EP" ss_type is specified, then filter taps will be moving
% alternately forward/backward by "n" units.  The reference should
% be the midpoint of the spectral lobes.    
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
% This method is also described in:
% (1) C.H. Cunningham, D.B. Vigneron, A.P. Chen, D. Xu, M. Lustig, D.A. Kelley,
% J.M. Pauly, Spectral?spatial excitation and refocusing for reduced volume
% mis- registration at 7 Tesla, in: Proceedings of the 14th Annual Meeting
% of ISMRM, Seattle, 2006, p. 72.
% (2) C.H. Cunningham, A.P. Chen, M. Lustig, J. Lupo, D. Xu, J. Kurhanewicz,
% R.E. Hurd, J.M. Pauly, S.J. Nelson, D.B. Vigneron, Pulse sequence for
% dynamic volumetric imaging of hyperpolarized metabolic products, J. Magn.
% Reson. 193 (1) (2008) 139?146.


% Error checking on inputs
    %
    if nargin < 7, 
	error('Usage: ss_spect_correct(b, bsf, N, n, f, ptype, ss_type, reg_factor, dbg)');
    end;
    
    if nargin < 8, 
	reg_factor = 0;
    end;
    
    if nargin < 9, 
	dbg = 0;
    end;
    
    switch ss_type, 
     case {'Flyback', 'EP'}
     otherwise
      error('ss_type must be ''Flyback'' or ''EP''');
    end;
    
    if length(bsf) ~= length(Noff), 
	error('bsf / offset vectors not same length');
    end;
    nfilt = length(bsf);

    % Calculate sampling positions of ref filter taps
    %
    N = length(b);
    t_ref = [0:N-1];
    
    % Get sampling in pass/stopbands
    % ensuring that mult_factor*N samples exist
    %
    mult_factor = 15;
    fdiff = diff(f);
    fsum = sum(fdiff(1:2:end));		% frequency range that is sampled
    df = fsum/(mult_factor*N);

    nband = length(f)/2;
    if strcmp(ss_type, 'Flyback'), 
	w = linspace(-pi, pi, 2*mult_factor*N);
    else
	w = [];
	for band = 1:nband,
	    nf = ceil((f(band*2)-f(band*2-1))/df) + 1;
	    df_act = (f(band*2)-f(band*2-1))/(nf-1);
	    wband = f(band*2-1) + [0:nf-1]*df_act;
	    w = [w pi*wband];
	end;
    end;

    % Plot w sampling
    %
    if (dbg >= 2),			% Verbose
	figure;
	subplot(2,1,1);
	plot(w/pi,ones(length(w)),'bx');
	hold on;
	for band = 1:nband, 
	    plot(f(band*2-1)*ones(1,2), [0 1], 'r');
	    plot(f(band*2)*ones(1,2), [0 1], 'r');
	end;
	title('Band Sampling');
	xlabel('Normalized Frequency');
	axis([min(w/pi) max(w/pi) -0.2 1.2]);

	[h_freq,freq] = freqz(b,1,[min(w/pi):0.005:max(w/pi)],2);
	
	subplot(2,1,2);
	plot(freq,abs(h_freq));
	title('Frequency Response');
	xlabel('Normalized Frequency');
	axis([min(w/pi) max(w/pi) -0.2 1.2]);
    end;
    

    % Calculate spectral filter corrections
    %
     if (dbg >= 2), 
	filt_fig = figure;
    end;
    
    % Get reference transform
    %
    Wref = exp(-i*kron(w', t_ref));
    Fref = Wref * b(:);

    for idx = 1:nfilt, 
	% Get actual sampling positions
	%
	switch (ss_type),
	 case 'Flyback', 
	  t_act = t_ref + Noff(idx)/Nper;
	 case 'EP', 
	  t_act = t_ref + (Noff(idx)/Nper * (-1).^[0:N-1]);
	end;
	
	if (dbg >= 2)
	    figure(filt_fig);
	    subplot(411);
	    stem([t_ref.' t_act.'], ones(length(t_ref),2));
	    legend('Reference', 'Actual');
	    title('Sampling Locations');
	end;

	switch (ss_type)
	 case 'Flyback'
	  % Get actual
	  %
	  Wact = exp(-i*kron(w', t_act));
	
	  % Get least-squares fit to filter 
	  %
	  type = 0;
	  switch(type)
	   case 0
	    % Least-squares solution with pseudo-inversion
	    % no regularization in this problem, typically not needed for Flyback designs
	    %
	    Wact_pinv = pinv(Wact);
	    bm(:,idx) = Wact_pinv * Fref;
	    %bm(:,idx) = inv(Wact'*Wact + 1e-4*eye(size(Wact,2)))*Wact'*Fref;
	   case 1
	    % Do constrained least-squares - best option, but slow
	    % Could be sped up using pdco
	    %
	    bm(:,idx) = lscon(Wact, Fref(:), 0, 1.2*max(abs(b)), b(:), 0);
	   case 2
	    % pdco option
	    %
	    c = 1;
	    pdco_options = pdcoSet;
	    d1 = 1e-6;
	    d2 = 1;
	    x0 = b(:);
	    y0 = zeros(size(Wact,1),1);
	    z0 = ones(size(Wact,2),1);
	    xsize = mean(abs(b));
	    zsize = 1;
	    bm(:,idx) = pdco(c, Wact, Fref(:), -1.2*max(abs(b)), 1.2*max(abs(b)), ...
			     d1, d2, pdco_options, x0, y0, z0, xsize,zsize,1);
	  end;
	  
	  if 0
	      figure;
	      Fref_fix = Wact * bm(:,idx);
	      plot(w/pi,abs(Fref_fix));
	      pause;
	  end;

	  % Since samples are still uniform, we can 
	  % just use SLR to get rf
	  %
	  if slr, 
	      rfm(:,idx) = b2rf(bsf(idx) * bm(:,idx));
	  else
	      rfm(:,idx) = 2*asin(abs(bsf(idx))) * conj(bm(:,idx));
	  end;
	 case 'EP'
	  % Get actual
	  %
	  Wact = exp(-i*kron(w', t_act));
	
	  % Get least-squares fit to filter 
	  %
	  %  Wact_pinv = pinv(Wact);
	  % bm(:,idx) = Wact_pinv * Fref;

      	  % Changed to regularized least-squares, works pretty well to keep peak B1 from
      	  % getting too large
	  bm(:,idx) = inv(Wact'*Wact + reg_factor*eye(size(Wact,2)))*Wact'*Fref;

	  % not very successful with constrained least-squares for keeping peak from getting large
      	  %   bm(:,idx) = lsqlin(Wact, Fref(:), [], [], [], [], zeros(size(b)), 2*max(abs(b))*ones(size(b)), b(:));
	  
	  % Can't do SLR yet
	  %
	  if slr, 
	      rfm(:,idx) = 2*asin(abs(bsf(idx)))*conj(bm(:,idx));
	  else
	      rfm(:,idx) = 2*asin(abs(bsf(idx))) * bm(:,idx);
	  end;
	end;
	
	if ( (dbg >= 2) && (rem(idx,10)==0) ), 
	    figure(filt_fig);
	    show_M = 0;
	    plot_db = 0;
	    if show_M,
		rf_ref = b2rf(bsf(round(end/2)) * b(:));
		subplot(412);
		hold off;
		plot(abs(rf_ref),'b-');
		hold on;
		plot(abs(rfm(:,idx)),'r--');

		% Fill in large time indices with RF
		%
		t_ref_x = Nper + t_ref*Nper;
		t_act_x = Nper + round(t_act*Nper);

		% Fill in large (mostly zero) arrays with RF
		%
		rf_ref_x = zeros(1,(N+2)*Nper);
		rf_ref_x(t_ref_x) = rf_ref;

		rf_act_x = zeros(1,(N+2)*Nper);
		rf_act_x(t_act_x) = rf_ref;

		rf_actfix_x = zeros(1,(N+2)*Nper);
		rf_actfix_x(t_act_x) = rfm(:,idx);
		
		% Now determine M 
		%
		g = pi * ones(length(rf_ref_x));
		Mref = ab2ex(abr(rf_ref_x, g, w/pi));
		Mact = ab2ex(abr(rf_act_x, g, w/pi));
		Mactfix = ab2ex(abr(rf_actfix_x, g, w/pi));
		subplot(413);
		if plot_db, 
		    hold off;
		    plot(w/pi,20*log10(abs(Mref)),'b');
		    hold on;
		    plot(w/pi,20*log10(abs(Mactfix)),'r--');
		    plot(w/pi,20*log10(abs(Mact)),'g--');
		    ylabel('DB Scale');
		else
		    hold off;
		    plot(w/pi,abs(Mref),'b');
		    hold on;
		    plot(w/pi,abs(Mactfix),'r--');
		    plot(w/pi,abs(Mact),'g--');
		    ylabel('Linear Scale');
		end;
		title('Magnitude Magnetization');
		
		subplot(414);
		hold off;
		plot(w/pi,angle(Mref),'b');
		hold on;
		plot(w/pi,angle(Mactfix),'r--');
		plot(w/pi,angle(Mact),'g--');
		title('Phase Magnetization');
	    else, 
		Fact = Wact * b(:);
		Fact_fix = Wact * bm(:,idx);

		subplot(4,1,2);
		hold off;
		plot(abs(b));
		hold on;
		plot(abs(bm(:,idx)), 'r--');
		title('Beta Polynomials');

		subplot(4,1,3);
		if plot_db, 
		    hold off;
		    plot(w/pi,20*log10(abs(Fref)),'b-');
		    hold on;
		    plot(w/pi, 20*log10(abs(Fact)), 'g--');
		    plot(w/pi, 20*log10(abs(Fact_fix)), 'r--');
		    ylabel('DB Scale');
		else
		    hold off;
		    plot(w/pi,abs(Fref),'b-');
		    hold on;
		    plot(w/pi, abs(Fact), 'g--');
		    plot(w/pi, abs(Fact_fix), 'r--');
		    ylabel('Linear Scale');
		end;
		
		title('Magnitude Response');

		subplot(4,1,4);
		hold off;
		plot(w/pi,angle(Fref),'b-');
		hold on;
		plot(w/pi,angle(Fact),'g--');
		plot(w/pi,angle(Fact_fix),'r--');
		title('Phase Response');
	    end;

	    fprintf(1,'Offset: %f   -- Hit any key to continue\n', Noff(idx));
	    pause;
	end;
    end;
    
