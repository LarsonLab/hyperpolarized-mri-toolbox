
function [h, status] = fir_pm_minpow(n, f, a, d, a_min, dbg)
% FIR_PM_MINPOW - FIR filter design using Parks-McLellan algs, 
%                 with attempt to minimize power     
%   
% Design n-tap linear-phase filter that meets multiband frequency 
% specification.  
%
% function [h, status] = fir_pm_minpow(n, f, a, d, a_min, dbg)
%    
% Inputs: 
%   n: number of taps returned
%   f: frequency bands
%   a: amplitude at band edges
%   d: ripple  in bands
%   a_min: if not present or [], chooses min(0, min(a-d))    
%   dbg: flag to turn on debugging statements/plots
%
% First designs filter with PM algorithm, then uses result as 
% a seed for the fir_linprog routine.  fir_linprog() includes
% a minimimization of the sum transition band response---which
% is equal to power when operating on an autocorrelation filter
% since the frequency response is the power spectrum (real and all    
% positive). When operating on the amplitude response, the 
% transition frequency response can both be negative and is 
% only proportional to the sqrt(power), so care should be taken
% when using it in this manner.
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
% $Header: /home/adam/cvsroot/src/ss/fir_pm_minpow.m,v 1.6 2012/02/01 00:41:22 peder Exp $
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
    oversamp = 16;

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
    wt = max(dn) ./ dn;
    lgrid = 31;				% Oversample, default 25
    try
	[h,d_opt,opt] = cfirpm(n-1,fn/pi,an,wt,{lgrid});
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
	resp_ok = check_response(fn/pi, an, dn, opt.fgrid, abs(opt.H));
    end;
    
    if (~resp_ok) 
	status = 'Failed';
	h = [];
	return;
    end;
    
    % Now call fir_linprog with designed filter as starting point.
    % The frequency response of the returned filter will be used
    % to refine our transition bands.
    if (dbg)
	fprintf(1,'Getting linear filter based on PM design\n');
    end;
    hlin = fir_linprog(n, f/pi, a, d, h(:), dbg);
    Hlin = freqz(hlin,1,w/pi,2);
    
    % Update transition bands
    %
    fn = [];
    an = [];
    dn = [];
    for tran = 1:ntran, 
	if tran == 1, 
	    f_l = min(w);		% This avoids sample at -pi,0 
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

	% Break up this transition region into bands
	% each with a max bound determined by the Hlin
	% filter response
	%
	% cfirpm seems to choke sometimes---I hypothesize
	% this is because bands are too close together, so
	% keep at least "nskip" points between neighbours
	%
	nskip = 1;
	ntran_pts = oversamp*2;
	ntran_region = max(1,round((length(idx_tran)-nskip)/ntran_pts));
	ntran_pts = (length(idx_tran)-nskip)/ntran_region;
	idx_region_st = nskip + 1 + round([0:ntran_region-1]*ntran_pts);
	idx_region_end = [idx_region_st(2:end) (length(idx_tran)-nskip+1)]-1;
	
	for reg = 1:length(idx_region_st)
	    % Get frequency indices corresponding to this region
	    %
	    idx_tran_reg = idx_tran(idx_region_st(reg):idx_region_end(reg));
	    if length(idx_tran_reg) == 0, 
		f_tran = [];
		a_tran = [];
		d_tran = [];
	    else
		f_tran = [min(w(idx_tran_reg)) max(w(idx_tran_reg))];
		max_H = max(abs(Hlin(idx_tran_reg))); % Phase not removed
						      % from Hlin
		tol_H = 0.01;
		max_H = max_H + tol_H*max(a+d2); % Add some tolerance to bands
		max_H = max(min(a+d2), max_H); % No tighter than stopband
		min_H = a_min;
		amp_tran = (max_H + min_H)/2;
		ripple_tran = (max_H - min_H)/2;
		a_tran = [amp_tran amp_tran];
		d_tran = [ripple_tran];
	    end;
	    fn = [fn f_tran];
	    an = [an a_tran];
	    dn = [dn d_tran];
	end;
	
	if tran < ntran,
	    fn = [fn f(tran*2-1) f(tran*2)];
	    an = [an a(tran*2-1) a(tran*2)];
	    dn = [dn d(tran)];
	end;
    end;	
	
    % Plot new frequency specfication on top of linear filter response
    %
    if (dbg >= 2) 
	figure;
	plot(w/pi, abs(Hlin));
	hold on;
	plot(opt.fgrid, abs(opt.H),'r');
	nband = length(fn)/2;
	for band = 1:nband, 
	    idx = [band*2-1:band*2];
	    plot(fn(idx)/pi, an(idx)+dn(band), 'g--');
	    plot(fn(idx)/pi, an(idx)-dn(band), 'g--');
	end;
	nband = length(f)/2;
	for band = 1:nband, 
	    idx = [band*2-1:band*2];
	    plot(f(idx)/pi, a(idx)+d(band), 'k--');
	    plot(f(idx)/pi, a(idx)-d(band), 'k--');
	end;
	fprintf(1,'Pausing...\n');
	pause;
    end;
    
    % Determine error weights, then call firpm
    %
    wt = max(dn) ./ dn;
    lgrid = 31;				% Oversample, default 25
    try
	[h,d_opt,opt] = cfirpm(n-1,fn/pi,an,wt,{lgrid});
    catch
	h = [];
	lsterr = lasterror;
	fprintf(1,'Error caught in cfirpm: \n');
	fprintf(1,'%s\n', lsterr.message);
    end;

    if (dbg >= 2) 
	plot(opt.fgrid, abs(opt.H),'m');
	fprintf(1,'Pausing...\n');
	pause;
    end;
    
    % Check frequency response at extremal frequencies
    % that they are within specified bands
    %
    resp_ok = 0;
    if ~isempty(h)
	resp_ok = check_response(fn/pi, an, dn, opt.fgrid, abs(opt.H));
    end;
    
    if (~resp_ok) 
	fprintf(1,'*** Failed to get min energy pulse ***\n');
	status = 'Failed';
	h = [];
	return;
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
