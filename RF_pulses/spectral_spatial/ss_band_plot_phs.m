function ss_band_plot_phs(f, a, d, fo, fs, f1, f2, sym)
% SS_BAND_PLOT - Plot filter specification
%   
% Inputs: 
% f - frequency band edges in normalized frequency [-1..1]
% a - amplitude of bands
% d - ripple weighting
% fo - offset absolute frequency
% fs - sampling frequency
% f1, f2 - upper lower absolute frequencies to plot
% sym - flag whether response is symmetric
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
% $Header: /home/adam/cvsroot/src/ss/ss_band_plot_phs.m,v 1.2 2013/08/15 03:34:50 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Exit if frequency vector is empty
    %
    if isempty(f), 
	fprintf(1,'No frequency bands to plot\n');
	return;
    end;
    
    if sym, 
	fmod = f;
	f = [-fmod(end:-1:1) fmod];
	a = [a(end:-1:1) a];
	d = [d(end:-1:1) d];
    end;

    % Plot unaliased frequency bands
    %
    nf = length(f)/2;
    figure;
    subplot(211);
    hold on;
    for idx = 1:nf, 
	plot(f((idx*2-1):(idx*2))*fs/2+fo, abs([a(idx) a(idx)]));
	plot(f((idx*2-1):(idx*2))*fs/2+fo, abs([a(idx) a(idx)]) + abs(d(idx)), '--');
	plot(f((idx*2-1):(idx*2))*fs/2+fo, abs([a(idx) a(idx)]) - abs(d(idx)), ...
	     '--');
    end;
    grid;
    title('Unaliased Filter Specification');
    xlabel('Frequency');
    ylabel('Amplitude Response');
    
    % Plot absolute frequency response
    %
    subplot(212);
    hold on;

    % Find what repeating responses need to be considered
    %
    low_mult = ceil(((fo-f1)-fs/2) / fs);
    hi_mult = ceil(((f2-fo)-fs/2) / fs);
    for mult = -low_mult:hi_mult, 
	foff = fo + fs*mult;
	for idx = 1:nf, 
	    % Get repeating band edge
	    %
	    f_l = (f(idx*2-1) * fs/2) + foff;
	    f_u = (f(idx*2) * fs/2) + foff;
	    
	    % Plot if in bounds
	    %
	    f_l = min(max(f_l, f1), f2);
	    f_u = min(max(f_u, f1), f2);
	    if (f_l ~= f_u), 
		plot([f_l f_u], abs([a(idx) a(idx)]));
		plot([f_l f_u], abs([a(idx) a(idx)]) + abs(d(idx)), '--');
		plot([f_l f_u], abs([a(idx) a(idx)]) - abs(d(idx)), '--');
	    end;
	end;
    end;
    grid;
    title(sprintf('Filter Specification - Fs: %6.1f Hz', fs));
    xlabel('Frequency');
    ylabel('Amplitude Response');

    
    
    
