% Based on 
% "FIR Filter Design via Spectral Factorization and Convex Optimization"
% by S.-P. Wu, S. Boyd, and L. Vandenberghe
%
% Determines minimum-order linear-phase filter that meets 
% amplitude/ripple specifications.
%
% function [h, status] = fir_min_order(n, f, a, d, dbg)
%
% Inputs: --- similar to cfirpm
%   n: max number of taps to try
%   f: frequency bands
%   a: amplitude at band edges
%   d: ripple  in bands
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
% $Header: /home/adam/cvsroot/src/ss/fir_min_order_qprog_phs.m,v 1.2 2012/06/27 16:49:53 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Bisection search calling fir_qprog to determine minimum order
% filter that meets linear phase requirements. By default it will check both
% odd and even length filters.  If a non-zero point exists at +/- 1, it will
% only check odd filters.  Only odd or even filters can be specified by the
% parameter odd_even.
%
% Inputs:
%   n - number of taps
%   f - frequency bands
%   a - amplitudes at band edges
%   even_odd: 1 - odd only
%             2 - even only
%             x - either odd or even
%   dbg - level of debug info to print
%

function [h, status] = fir_min_order_qprog_phs(n, f, a, d, even_odd, dbg)
    if nargin < 4, 
	error(['Usage: function [h, status] = fir_min_order(n, f, a, d, even_odd,' ...
	       ' dbg)']);
    end;
    
    status = 'Failed';
    h = [];
    
    if nargin < 5 || isempty(even_odd), 
	even_odd = 0;
    end;
    switch (even_odd)
     case {1, 2}
     otherwise
      even_odd = 0;
    end;
    
    if nargin < 6, 
	dbg = 0;
    end;
    
    % Initialize best odd/even filters
    %
    hbest_odd = [];
    hbest_even = [];
    
    % Get max odd/even num taps
    %
    n_odd_max = 2*floor((n-1)/2)+1;
    n_even_max = 2*floor(n/2);

    % Test odd filters first 
    %
    if even_odd ~= 2, 
	n_bot = 1;
	n_top = (n_odd_max+1)/2;
	n_cur = n_top;
	if (dbg)
	    fprintf(1, 'Testing odd length filters...\n');
	end;
	while (n_top - n_bot > 1), 
	    n_tap = (n_cur * 2 - 1);
	    if (dbg)
		fprintf(1, '%4d taps: ...', n_tap)
	    end;
	    [h, status] = fir_qprog_phs(n_tap, f, a, d, hbest_odd, dbg);
	    if strcmp(status, 'Solved')
		% feasible
		hbest_odd = h;
		if dbg, 
		    fprintf(1,'Feasible\n');
		end;
		n_top = n_cur;
		if n_top == n_bot+1, 
		    n_cur = n_bot;
		else
		    n_cur = ceil((n_top + n_bot)/2);
		end;
	    else
		if dbg, 
		    fprintf(1,'Infeasible\n');
		end;
		n_bot = n_cur;
		n_cur = ceil((n_bot+n_top)/2);
	    end;
	end;
    end

    % Test even filters now
    %
    if even_odd ~= 1, 
	n_bot = 1;
	if isempty(hbest_odd), 
	    n_top = n_even_max/2;
	    n_cur = n_top;
	else
	    n_top = min(n_even_max/2, (length(hbest_odd)+1)/2);
	    n_cur = n_top;
	end;
	if (dbg)
	    fprintf(1, 'Testing even length filters...\n');
	end;
	while (n_top - n_bot > 1), 
	    n_tap = n_cur * 2;
	    if (dbg)
		fprintf(1, '%4d taps: ...', n_tap)
	    end;
	    [h, status] = fir_qprog_phs(n_tap, f, a, d, hbest_even, dbg);
	    if strcmp(status, 'Solved')
		% feasible
		hbest_even = h;
		if dbg, 
		    fprintf(1,'Feasible\n');
		end;
		n_top = n_cur;
		if n_top == n_bot+1, 
		    n_cur = n_bot;
		else
		    n_cur = ceil((n_top + n_bot)/2);
		end;
	    else
		if dbg, 
		    fprintf(1,'Infeasible\n');
		end;
		n_bot = n_cur;
		n_cur = ceil((n_bot+n_top)/2);
	    end;
	end;
    end
    
    if isempty(hbest_odd) && isempty(hbest_even), 
	status = 'Failed';
	h = [];
	if dbg, 
	    fprintf(1,'\nFailed to achieve specs\n');
	end;
    else
	status = 'Solved';
	if (isempty(hbest_odd)) 
	    h = hbest_even;
	elseif (isempty(hbest_even))
	    h = hbest_odd;
	elseif (length(hbest_odd) < length(hbest_even))
	    h = hbest_odd;
	else
	    h = hbest_even;
	end;
	
        %  dbg = 2;
        if (dbg >= 2) 
            filt_fig = figure;

            m = 512;
            H = fftf(h, m);
            freq = [-m/2:m/2-1]/m*2;

            % Correct H by half-sample offset if even number of taps
            %
            n_tap = length(h);
            if bitget(n_tap,1) == 0, 
                H = H .* exp(-i*pi*freq(:)*0.5);
            end;

            figure(filt_fig);
            clf;
            hold on;
            plot_spec_phs(f,a,d,freq,H);
            title('Frequency Response');
            xlabel('Normalized Frequency');
            fprintf(1,'Pausing...');
            pause;
            fprintf(1,'\r           \r');
        end;

        if dbg, 
	    fprintf(1,'\nOptimum number of filter taps is: %d.\n',length(h));
	end;
    end;

