
function [h, status] = fir_pm(n, f, a, d, a_min, dbg)
% FIR_LINPROG - FIR filter design using Parks-McLellan algs
%   
% Design n-tap linear-phase filter that meets multiband frequency 
% specification.  
%
% function [h, status] = fir_pm(n, f, a, d, dbg)
%    
% Inputs: 
%   n: number of taps returned
%   f: frequency bands
%   a: amplitude at band edges
%   d: ripple  in bands
%   a_min: if not present, chooses min(0, min(a-d))    
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
% $Header: /home/adam/cvsroot/src/ss/fir_pm.m,v 1.7 2012/02/01 00:41:22 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    % Default value for a_min
    %
    d2 = [d(:).'; d(:).'];
    d2 = d2(:).';
    if (nargin < 5) || isempty(a_min), 
	a_min = min(0,min(a-d2));
    end;

    % Default value for dbg
    %
    if nargin < 6, 
	dbg = 0;
    end;
    
    % Determine if real or complex coefficients
    %
    f = f * pi;				% Scale to +/- pi
    if min(f) < 0, 
	real_filter = 0;
    else
	real_filter = 1;
    end;
    
    % Determine if filter has odd or even number of 
    % taps
    %
    if (bitget(n,1) == 1) 
	odd_filter = 1;
    else
	odd_filter = 0;
    end;
    
    % If the frequency specification has a non-zero point
    % at +/- 1, then the order must be even. A warning is 
    % printed and a failure returned if this is the case.
    %
    if (~odd_filter)
	idx = find(abs(f) ~= 0);
	if find(a(idx) ~= 0)
	    warning('n odd and frequency spec non-zero at fs/2');

	    status = 'Failed';
	    h = [];
	    return;
	end;
    end;
    
    % Oversampling on frequency to determine transition bands
    %
    oversamp = 8;

    % Get first pass on w
    %
    if real_filter, 
	m = oversamp * n;
	w = linspace(0,pi,m);
    else
	m = 2 * oversamp * n;
	w = linspace(-pi,pi,m);
    end;

    % Find bounds on transition regions and convert to amp/ripple
    %
    ub_tran = max(a + d2);
    lb_tran = a_min;			% Set to min amplitude spec
    amp_tran = (ub_tran + lb_tran)/2;
    ripple_tran = (ub_tran - lb_tran)/2;
    
    % Find indices of transition bands, build up new frequency spec
    %
    nband = length(f)/2;
    ntran = nband+1;
    fn = [];
    an = [];
    dn = [];
    for tran = 1:ntran, 
	if tran == 1, 
	    f_l = min(w);		% This avoids sample at -pi 
	    rband = tran;
	    f_r = f(rband*2-1);
	elseif tran == ntran, 
	    lband = tran-1;
	    f_l = f(lband*2);
	    f_r = pi;			% This avoids sample at pi
	else
	    lband = tran-1;
	    f_l = f(lband*2);
	    rband = tran;
	    f_r = f(rband*2-1);
	end;
	idx_tran = find((w > f_l) & (w < f_r));
	% cfirpm seems to choke sometimes---I hypothesize
	% this is because the transition edges are too 
	% close to the actual passbands, so don't take 
	% the immediately adjacent points
	%
	nskip = 1;
	if length(idx_tran) <= 1+2*nskip,
	    f_tran = [];
	    a_tran = [];
	    d_tran = [];
	else
	    idx_tran = idx_tran(1+nskip:end-nskip);
	    f_tran = [min(w(idx_tran)) max(w(idx_tran))];
	    a_tran = [amp_tran amp_tran];
	    d_tran = [ripple_tran];
	end;
	fn = [fn f_tran];
	an = [an a_tran];
	dn = [dn d_tran];
	if tran < ntran,
	    fn = [fn f(tran*2-1) f(tran*2)];
	    an = [an a(tran*2-1) a(tran*2)];
	    dn = [dn d(tran)];
	end;
    end;	

    % Determine error weights, then call firpm
    %
    w = max(dn) ./ dn;
    lgrid = 31;				% Oversample, default 25
    if 0
	% firpm has some instability but cfirpm seems ok...
	%
    if real_filter, 
	try
	    [h,d_opt,opt] = firpm(n-1,fn/pi,an,w,{lgrid});
	catch
	    h = [];
	end;
    else
	[h,d_opt,opt] = cfirpm(n-1,fn/pi,an,w,{lgrid});
    end;
    end;

%    [h,d_opt,opt] = cfirpm(n-1,fn/pi,an,w,{lgrid},'skip_stage2');
    try
	[h,d_opt,opt] = cfirpm(n-1,fn/pi,an,w,{lgrid});
    catch
	h = [];
	lsterr = lasterror;
	fprintf(1,'Error caught in cfirpm: \n');
	fprintf(1,'%s\n', lsterr.message);
    end;

    % Check frequency response at extremal frequencies
    % that are within specified bands
    %
    resp_ok = 0;
    if ~isempty(h)
	resp_ok = check_response(f/pi, a, d, opt.fgrid, abs(opt.H));
    end;
    
    if (~resp_ok) 
	status = 'Failed';
	if dbg>=2, 
	    plot_response(opt.fgrid, opt.H, fn/pi, an, dn);
	    title('Filter Response');
	    pause(1);
	end;
        h = [];
    else
	if dbg>=2, 
	    plot_response(opt.fgrid, opt.H, fn/pi, an, dn);
	    title('Filter Response');
	    pause(1);
	end;

	h = h(:);
	status = 'Solved';
    end;
return;
    
    
function status = check_response(f,a,d,ftest,htest)
% CHECK_RESPONSE - Check magnitude response to see if it meets specs
%   
    nband = length(f)/2;
    status = 1;
    for band = 1:nband, 
	idx = find((ftest >= f(band*2-1)) & (ftest <= f(band*2)));
	if isempty(idx)
	    break;
	end;

	f_off = ftest(idx) - f(band*2-1);
	a_test = a(band*2-1) + ...
		 (a(band*2)-a(band*2-1)) * f_off/(f(band*2)-f(band*2-1));

	a_hi = a_test + d(band);
	a_lo = a_test - d(band);

	if (find((htest(idx) > a_hi) | (htest(idx) < a_lo)))
	    status = 0;			% Fails in at least one sample
	    return;
	end;
    end;
    
    return;
    
function plot_response (freq,h,f,a,d)
% plot_response - Plot frequency specification and actual response
%   
    figure;
    hold on;
    nband = length(f)/2;
    for band = 1:nband, 
	idx = [band*2-1:band*2];
	plot(f(idx), a(idx)+d(band), 'k--');
	if max(a(idx)-d(band)) > 0, 
	    plot(f(idx), max(0,a(idx)-d(band)), 'k--');
	end;
    end;
    plot(freq, real(h));
    plot(freq, imag(h),'b--');
    xlabel('Frequency');
    ylabel('Filter Response');

    return;
