% function [h, fn, an, dn] = fir_expand(n, f, a, d, a_min, dbg)
%
% Inputs: 
%   n: number of taps
%   f: frequency bands
%   a: amplitude at band edges
%   d: ripple in bands
%   a_min: minimum amplitude to allow, default if empty is min(0,min(a-d));
%   dbg: flag to turn on debugging statements/plots
%
% Outputs
%   h: new filter 
%   fn: new frequency bands
%   an: new amplitudes
%   dn: new ripple
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
% $Header: /home/adam/cvsroot/src/ss/fir_expand.m,v 1.5 2012/02/01 00:41:22 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [hn,fn,an,dn,status] = fir_expand(n,f,a,d,a_min,dbg)

    % Check input parameters 
    %
    if nargin < 4, 
	error(['Usage: function [h,fn,an,dn,status] = fir_expand(n,f,a,' ...
	       'd,[dbg])']);
    end;
    
    if nargin < 5, 
	a_min = [];
    end;

    if nargin < 6, 
	dbg = 0;
    end;
    
    if (find(a(1:2:end)-a(2:2:end)))
	error('Only works for flat bands now');
    end;
    
    % Confirm baseline design works
    %
    fn = f;
    an = a;
    dn = d;
    [hn,status] = fir_pm(n,f,a,d,a_min,dbg);
    if strcmp(status,'Failed')
%	warning('Baseline design infeasible.');
	return;
    end;
    
    % Now get number of transition regions, and type of filter
    %
    nband = length(f)/2;
    if min(f) < 0, 
	real_filter = 0;
    else
	real_filter = 1;
    end;
    
    if bitget(n,1) == 1, 
	odd_filter = 1;
    else
	odd_filter = 0;
    end;
    
    % Get frequency grid
    %
    oversamp = 7;
    if real_filter
	m = oversamp * n;
	fg = [0:m-1]/(m-1);
    else
	m = 2 * oversamp * n;
	fg = [-m/2:m/2-1]/(m/2);
    end;
    
    % Keep expanding each filter band up to point filter fails
    % Choose lower-amplitude band to expand if possible, else expand
    % bands equally
    %
    fbest = f;
    if real_filter, 
	ntran = nband+1;
    else
	ntran = nband;
    end;
    if dbg >= 2, 
	ftry_fig = figure;
    end
    
    for tran = 1:ntran, 
	if dbg >= 1, 
	    fprintf(1, 'Optimizing transition %d of %d\n', tran, ntran);
	end;
	
	% Get indices of left/right bands
	%
	if real_filter 
	    if (tran == 1) 
		lband = 0;		% Real-filter, no band left of band 1
		rband = 1;
	    elseif (tran == nband+1)
		lband = nband;
		rband = 0;
	    else 
		lband = tran-1;
		rband = tran;
	    end;
	else
	    lband = mod(tran-2,nband)+1;
	    rband = mod(tran-1,nband)+1;
	end;

	% Get starting fgrid indices corresponding to 
	% band edges.
	%
	if lband ~= 0,
	    lidx = find((fg <= f(lband*2)), 1, 'last');
	end;
	if rband ~= 0, 
	    ridx = find((fg >= f(rband*2-1)), 1, 'first');
	end;
	
	% Get max number of grid points the band edges 
	% can move, as well as the step size for each 
	% band
	%
	if rband == 0, 
	    nlmax = m - lidx;
	    nrmax = 0;
	elseif lband == 0, 
	    nlmax = 0;
	    nrmax = ridx - 1;
	else
	    nlmax = mod(ridx-lidx, m);	
	    nrmax = nlmax;
	end;
	ntotal = max(nlmax, nrmax);
	
	% Find order of band edges to expand
	% -- preference is to expand stopbands first, then passbands
	%    in order to minimize energy in filter

	% Only add bands that are stopbands
	%
	band_order = [];
	if ((lband ~= 0) && (rband~=0))
	    if a(lband*2) == min(a) 
		band_order = [0];
	    end;
	    if (a(rband*2) == min(a))
		band_order = [band_order 1];
	    end;
	elseif (lband ~= 0)
	    if a(lband*2) == min(a) 
		band_order = [band_order 0];
	    end;
	elseif (rband ~= 0)
	    if (a(rband*2) == min(a))
		band_order = [band_order 1];
	    end;
	end;

	if 0
	band_order = [0 1];		% 0 = left, 1 = right
	if ((lband ~= 0) && (rband ~= 0) && ...
	    (a(lband*2) > a(rband*2-1))), % swap order
	    band_order = [1 0];
	end;
	end;
	
	    
	% Bisection search on how many points to add to band
	% edges. nlbot and nrbot are known working solutions.
	%
	nlbot = 0; nltop = nlmax; 
	nrbot = 0; nrtop = nrmax; 
	for working_band = band_order, 
	    if working_band == 0,	% lband
		if dbg >= 2,
		    fprintf(1, 'Expanding left edge...\n');
		end;
		% Update nltop for any right edge expansion
		%
		if lband ~= 0 
		    nltop = min(nltop, ntotal-nrbot);
		end
		while ((nltop - nlbot) > 0)

		    if dbg >= 2, 
			fprintf(1,'nlbot: %d  nltop: %d -- nrbot: %d  nrtop: %d -- ntotal: %d\n', ...
				nlbot, nltop, nrbot, nrtop, ntotal);
		    end;
		    % Advance left band and calculate feasibility
                    % 
			nltry = ceil((nltop + nlbot)/2);
			lidxtry = mod(lidx + nltry -1, m) + 1;
			ftry_cont = fbest;
			ftry_cont(lband*2) = fg(lidxtry);

			% Create new freq spec from ftry_cont that removes
			% any possible wraps (only occurs for complex filters)
			%
			ftry = ftry_cont;
			atry = a;
			dtry = d;
			if ftry(1) > ftry(2), 
			    if (ftry(1) == 1), 
				ftry(1) = -1;
			    else
				ftry = [-1 ftry(2:end) ftry(1) 1];
				atry = [a a(1) a(2)];
				dtry = [d d(1)];
			    end;
			elseif ftry(end-1) > ftry(end), 
			    if (ftry(end) == -1)
				ftry(end) = 1;
			    else
				ftry = [-1 ftry(end) ftry(1:end-1) 1];
				atry = [a(end) a(end) a];
				dtry = [d(end) d];
			    end;
			end;
			
			if dbg>=2, 
			    clf;
			    hold on;
			    nb = length(ftry)/2;
			    for bnd = 1:nb, 
				plot([ftry(bnd*2-1) ftry(bnd*2-1) ftry(bnd*2) ftry(bnd*2)], ...
				     log10(1+[0 (atry(bnd*2-1)+dtry(bnd)) (atry(bnd*2)+dtry(bnd)) 0]), 'r');
			    end;
			    for bnd = 1:nband, 
				plot([f(bnd*2-1) f(bnd*2-1) f(bnd*2) f(bnd*2)], ...
				     log10(1+[0 (a(bnd*2-1)+d(bnd)) (a(bnd*2)+d(bnd)) 0]), 'b');
			    end;
			    drawnow;
			    fprintf(1,'<Pausing>');
			    pause;
			    fprintf(1,'\r         \r');
			end;

			if any(diff(ftry) <= 0)
			    status = 'Failed';
			else
			    [h,status] = fir_pm(n,ftry,atry,dtry,a_min,dbg);
			end;
			if strcmp(status, 'Solved')
			    if dbg >= 3, 
				fprintf(1,'nltry: %d Solved\n',nltry);
			    end;
			    % Update nlbot
			    %
			    nlbot = nltry;

			    % Update current best results
			    %
			    fbest = ftry_cont;

			    % Update output
			    %
			    hn = h;
			    fn = ftry;
			    an = atry;
			    dn = dtry;
			else
			    if dbg >= 3,
				fprintf(1,'nltry: %d Failed\n',nltry);
			    end;
			    % Update nltop
			    %
			    if nltop == nlbot+1
				nltop = nlbot;
			    else
				nltop = nltry;
			    end;
			end;
		end; % while
	    else % working band is rband
		if dbg >= 2,
		    fprintf(1, 'Expanding right edge...\n');
		end;
		% Update nrtop for any left edge expansion
		%
		if rband ~= 0 
		    nrtop = min(nrtop, ntotal-nlbot);
		end
		while ((nrtop - nrbot) > 0)
		    if dbg >= 2, 
			fprintf(1,'nlbot: %d  nltop: %d -- nrbot: %d  nrtop: %d -- ntotal: %d\n', ...
				nlbot, nltop, nrbot, nrtop, ntotal);
		    end;
		    % Advance right band and calculate feasibility
		    %
			nrtry = ceil((nrtop + nrbot)/2);
			ridxtry = mod(ridx - nrtry -1, m) + 1;
			ftry_cont = fbest;
			ftry_cont(rband*2-1) = fg(ridxtry);

			% Create new freq spec from ftry_cont that removes
			% any possible wraps (only occurs for complex filters)
			%
			ftry = ftry_cont;
			atry = a;
			dtry = d;
			if ftry(1) > ftry(2), 
			    if (ftry(1) == 1)
				ftry(1) = -1;
			    else
				ftry = [-1 ftry(2:end) ftry(1) 1];
				atry = [a a(1) a(2)];
				dtry = [d d(1)];
			    end;
			elseif ftry(end-1) > ftry(end), 
			    if (ftry(end) == -1)
				ftry(end) = 1;
			    else
				ftry = [-1 ftry(end) ftry(1:end-1) 1];
				atry = [a(end) a(end) a];
				dtry = [d(end) d];
			    end;
			end;
			
			if dbg>=2, 
			    clf;
			    hold on;
			    nb = length(ftry)/2;
			    for bnd = 1:nb, 
				plot([ftry(bnd*2-1) ftry(bnd*2-1) ftry(bnd*2) ftry(bnd*2)], ...
				     log10(1+[0 (atry(bnd*2-1)+dtry(bnd)) (atry(bnd*2)+dtry(bnd)) 0]), 'r');
			    end;
			    for bnd = 1:nband, 
				plot([f(bnd*2-1) f(bnd*2-1) f(bnd*2) f(bnd*2)], ...
				     log10(1+[0 (a(bnd*2-1)+d(bnd)) (a(bnd*2)+d(bnd)) 0]), 'b');
			    end;
			    drawnow;
			    fprintf(1,'<Pausing>');
			    pause;
			    fprintf(1,'\r         \r');
			end;

			if any(diff(ftry) <= 0)
			    status = 'Failed';
			else
			    [h,status] = fir_pm(n,ftry,atry,dtry,a_min,dbg);
			end;
			if strcmp(status, 'Solved')
			    if dbg >= 3,
				fprintf(1,'nrtry: %d Solved\n',nrtry);
			    end;

			    % Update nrbot
			    %
			    nrbot = nrtry;

			    % Update current best results
			    %
			    fbest = ftry_cont;

			    % Update output
			    %
			    hn = h;
			    fn = ftry;
			    an = atry;
			    dn = dtry;
			else
			    if dbg >= 3,
				fprintf(1,'nrtry: %d Failed\n',nrtry);
			    end;
			    % Update nrtop
			    %
			    if nrtop == nrbot+1
				nrtop = nrbot;
			    else
				nrtop = nrtry;
			    end;
			end;
		end; % while
	    end; % if working band
	end; % for working band
	if dbg >= 2, 
	    fprintf(1,'nlbot: %d  nltop: %d -- nrbot: %d  nrtop: %d -- ntotal: %d\n', ...
		    nlbot, nltop, nrbot, nrtop, ntotal);
	end;
    end % for tran
    
    status = 'Solved';

    
