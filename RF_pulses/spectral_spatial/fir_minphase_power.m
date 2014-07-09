% Determines minimum-order minimum-phase filter that has a magnitude
% response that meets the magnitude amplitude/ripple specifications.  
% It will only return filters with an odd number of taps.
%
% function [h, status] = fir_minphase_power(n, f, a, d, use_max, dbg)
%
% Inputs: --- similar to cfirpm
%   n: max number of taps to try
%   f: frequency bands
%   a: amplitude at band edges
%   d: ripple  in bands
%   use_max: don't search for min order, use n taps
%   dbg: flag to turn on debugging statements/plots
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
% $Header: /home/adam/cvsroot/src/ss/fir_minphase_power.m,v 1.7 2012/06/19 16:53:51 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calls fir_min_order to determine minimum-length odd filter that 
% meets magnitude-squared response, then performs spectral 
% factorization
%

function [hn,status] = fir_minphase_power(n, f, a, d, use_max, dbg)
    if nargin < 4, 
	error('Usage: function [h, status] = fir_minphase_power(n, f, a, d, use_max,dbg)');
    end;
    
    if nargin < 5, 
	use_max = 0;
    end;

    if nargin < 6, 
	dbg = 0;
    end;
    
    hn = [];
    status = 'Failed';

    % Get magnitude-squared spec
    %
    d2 = [d(:).'; d(:).'];
    d2 = d2(:).';
    mxspec = (a+d2).^2;
    mnspec = max(0,(a-d2)).^2;

    % Get "zero" threshold of 2% of lowest spec
    %
    ztol_perc = 2;
    ztol = min(ztol_perc/100 * (mxspec-mnspec));
    
    % Offset magnitude-squared spec by ztol, get n
    %
    mnspec = max(ztol,mnspec);
    a_sqr = (mxspec + mnspec)/2;
    d2_sqr = (mxspec-mnspec)/2;
    d_sqr = d2_sqr(1:2:end);
    n_max = 2*n - 1;
    
    % Get minimum-order linear-phase filter that meets
    % magnitude response
    %
    odd_only = 1;
    if ~use_max, 
        [r, status] = fir_min_order(n_max, f, a_sqr, d_sqr, odd_only, ztol, dbg);
    else
        [r, status] = fir_pm(n_max, f, a_sqr, d_sqr, ztol, dbg);
    end;
    if strcmp(status, 'Failed')
        fprintf(1,'Failed to get filter\n');
        return;
    end

    Rok = 0;
    while ~Rok, 
	% Get min power filter
	%
	[rn,status] = fir_pm_minpow(length(r), f, a_sqr, d_sqr, ztol, dbg);
	
	oversamp = 15;
	m = 2 * oversamp * length(rn);
	m2 = 2^ceil(log2(m));
	R = real(fftf(rn,m2));

	% Check magnitude response to make sure that it is everywhere positive
	%
	if dbg,
	    freq = [-m2/2:m2/2-1]/m2*2;
	    Ro = real(fftf(r,m2));	% Linear-phase must have real autocorrelation
	    figure;
	    plot(freq(:),real(Ro(:)));
	    hold on;
	    plot(freq(:),real(R(:)),'r');
	    plot_spec(f,a_sqr,d_sqr,'k');
	    title('Squared Frequency Response of Autocorrelation Fcns');
	    xlabel('Normalized Frequency');
	end;
	
	% Now use spectral factorization to get return filter
	% --- first offset to make sure it's positive
	if min(R) < 0, 
	    min_stop = min(a_sqr + d2_sqr);
	    Rtol_perc = 0.1;		% 10% stopband tolerance
	    Rtol = Rtol_perc * min_stop;		% Rule of thumb
	    if min(R) < -Rtol, 
		fprintf(1, 'Autocorrelation has negative value\n');
		fprintf(1, '  Tol (%d%% stopband): %e  Actual: %e\n', ...
			round(Rtol_perc*100), Rtol, -min(R));
		
		% Test spectral factorization
		%
		rn = rn + Rtol;
		hn = spectral_fact(rn);
		hn = conj(hn(end:-1:1));
		
		% Get squared frequency response and check against specs
		% + Rtol
		%
		H = (fftf(hn, m2));
		freq = [-m2/2:m2/2-1]/m2*2;
		H2 = abs(H).^2;		% Squared-mag response (H is Min-Phase)
		nband = length(f)/2;
		atol = 0.05;
		Rok = 1;
		for band = 1:nband, 
		    idx = find((freq >= f(band*2-1)) & (freq <= f(band*2)));
		    amax = (1+atol)*(a_sqr(band*2-1) + d_sqr(band) + Rtol);
		    amin = (1-atol)*(a_sqr(band*2-1) - d_sqr(band));
		    fail = find((H2(idx) > amax) | ...
				(H2(idx) < amin));
		    if fail, 
			fprintf(1, '  Spectral factorization doesn''t meet specs\n');
			fprintf(1, '  Increase number of taps to: %d\n', ...
				length(r)+2);
    	    if dbg, 
			fprintf(1,'<Pausing>');
			figure;
			plot_spec(f,a_sqr,d_sqr,'k');
			hold on;
			plot(freq(:),R(:),'r');
			plot(freq(:),H2(:));
			grid;
			
			title('Squared Frequency Response of Original and Factorized Filter');
			xlabel('Normalized Frequency');
			pause;
            end
			fprintf(1,'\r          \r');
			[r, stat] = fir_pm_minpow(length(r)+2, f, a_sqr, d_sqr, ztol, ...
					   dbg);
			Rok = 0;
			break;
		    end;
		end;
	    else
		    fprintf(1, 'Autocorrelation has negative value, but within tol\n');
		    fprintf(1, '  Tol (%d%% stopband): %e  Actual: %e\n', ...
			    round(Rtol_perc*100), Rtol, -min(R));
		rn = rn - min(R);
		Rok = 1;
	    end;
    elseif strcmp(status, 'Failed')
        fprintf(1,'\r          \r');
			[r, stat] = fir_pm_minpow(length(r)+2, f, a_sqr, d_sqr, ztol, ...
					   dbg);
		Rok = 0;

    else
	    if dbg, 
            fprintf(1,'Autocorrelation OK\n');
	    end;
	    Rok = 1;
	end;
    end;

    % Spectral factorize and time reverse
    %
    h = spectral_fact(r);
    h = conj(h(end:-1:1));

    hn = spectral_fact(rn);
    hn = conj(hn(end:-1:1));
    if dbg, 
	m = 512;
	freq = [-m/2:m/2-1]/m*2;
	H = abs(fftf(h,512));
	Hn = abs(fftf(hn,512));
	figure;
	plot(freq,H);
	hold on;
	plot(freq,Hn,'r');
	plot_spec(f,a,d,'k');
	title('Frequency Response Factorized Filters');
    end;
    
