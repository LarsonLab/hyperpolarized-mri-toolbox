
function [h, status] = fir_linprog(n, f, a, d, h0, dbg)
% FIR_LINPROG - FIR filter design using linear programming
%   
% Design n-tap linear-phase filter that meets multiband frequency 
% specification.  If h0 is provided, it is used as an initial 
% estimate of the filter taps. 
%
% Also constrain H(f) > min(0, min(a(1:2:end)-d), min(a(2:2:end)-d))
% Also minimizes sum(H(f)) in transition bands.
%
% function [h, status] = fir_linprog(n, f, a, d, h0, dbg)
%    
% Inputs: --- similar to cfirpm
%   n: number of taps returned
%   f: frequency bands (-1->1)
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
% $Header: /home/adam/cvsroot/src/ss/fir_linprog.m,v 1.11 2013/08/15 16:01:59 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
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
	idx = find(abs(f) == pi);
	if find(a(idx) == 1)
	    warning('n odd and frequency spec 1 at fs/2');

	    status = 'Failed';
	    h = [];
	    return;
	end;
    end;
    
    % Determine number of optimization parameters
    %
    nhalf = ceil(n/2);		% number of taps in half-side of
				% filter
    nx = nhalf;
    if ~real_filter, 
	if odd_filter, 
	    nx = 2*nhalf-1;
	else
	    nx = 2*nhalf;
	end;
    end;

    % Create optimization arrays
    %
    oversamp = 15;
    undersamp_tran = 1;		% Undersampling factor for transition
                                % regions
    % Get first pass on w
    %
    if real_filter, 
	m = oversamp * n;
	w = linspace(0,pi,m);
    else
	m = 2 * oversamp * n;
	w = linspace(-pi,pi,m);
    end;
    
    % Add explicit samples to w at the edge of each specified band
    %
    w = sort([w f]);

    % Find indices to passbands/stopbands, and fill in upper/lower bounds
    %
    idx_band = []; U_band = []; L_band = [];
    nband = length(f)/2;
    for band = 1:nband, 
	idx = find( (w >= f(band*2-1)) & (w <= f(band*2)) );
	% Get amplitude from linear interpolation on band
	%
	idx_band = [idx_band idx];
	if (f(band*2-1) == f(band*2))
	    amp = a(band*2-1);
	else
	    amp = a(band*2-1) + (a(band*2)-a(band*2-1)) * ...
		  ((w(idx) - f(band*2-1))/(f(band*2)-f(band*2-1)));
	end;
	U_band = [U_band (amp + d(band))];
	L_band = [L_band (amp - d(band))];
    end; 

    % Get transition indices
    %
    idx_tmp = ones(1,length(w));
    idx_tmp(idx_band) = 0;
    idx_tran = find(idx_tmp == 1);

    % Get average representation of response
    %
    lb_resp = size(w);
    lb_resp(idx_band) = (U_band + L_band)/2;
    lb_resp(idx_tran) = (max(U_band) + min(L_band))/2;
    if real_filter, 
	lb_resp = [lb_resp(end:-1:1) lb_resp(2:end-1)];
    end;
    if dbg >= 3, 
	if real_filter, 
	    wplot = [-w(end:-1:1) w(2:end-1)];	
	else
	    wplot = w;
	end;
	figure;
	plot(wplot,lb_resp);
    end;

    % Fill in x0 with whatever taps are available from h0
    % - In updated linprog() this may not be necessary
    if nargin < 5 || isempty(h0),
	h0 = [];
    end;
    x0 = fill_opt_param(h0,nx,real_filter,odd_filter,lb_resp,dbg);
    
    % Decimate w in transition regions
    %
    idx_tran = idx_tran(1:undersamp_tran:end);

    % Add transition band limits to be between the + max 
    % specification on each band and min of (0,min(L_band))
    %
    if ~isempty(idx_tran)
	U_amp_tran = max(U_band);
	U_tran = U_amp_tran*ones(1,length(idx_tran));
	L_amp_tran = min(0, min(L_band));
	L_tran = L_amp_tran*ones(1,length(idx_tran));
    else
	U_tran = [];
	L_tran = [];
    end;
    
    % Update w, idx_band
    %
    wband = w(idx_band);
    idx_band = [1:length(wband)];
    wtran = w(idx_tran);
    idx_tran = [1:length(wtran)] + length(wband);
    w = [wband(:).' wtran(:).'];
    m = size(w,2);

    if dbg >= 3, 
	figure;
	plot(w(idx_band),U_band,'*');
	hold on;
	plot(w(idx_band),L_band,'o');
	plot(w(idx_tran),U_tran,'r*');
	plot(w(idx_tran),L_tran,'ro');
	pause;
    end;

    if real_filter
	% create optimization matrices
	% A is the matrix used to compute the power spectrum
	% A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos(n*w)]
	if (odd_filter)
	    Acos = [ones(m,1) 2*cos(kron(w',[1:nhalf-1]))];
	else
	    Acos = [2*cos(kron(w',[0:nhalf-1]+0.5))];
	end;
	Asin = [];
    else
	if (odd_filter) 
	    Acos = [ones(m,1) 2*cos(kron(w',[1:nhalf-1]))];
	    Asin = [2*sin(kron(w',[1:nhalf-1]))]; 
	else
	    Acos = [2*cos(kron(w',[0:nhalf-1]+0.5))];
	    Asin = [2*sin(kron(w',[0:nhalf-1]+0.5))]; 
	end;
    end;

    % Get subset of A matrix for current order
    %
    A = [Acos Asin];
    
    % Build matrix for upper bound constraints
    %
    A_U = [A(idx_band,:); A(idx_tran,:)];
    U_b = [U_band U_tran];
    
    % Build matrices for lower bound constraints
    %
    A_L = [A(idx_band, :); A(idx_tran,:)];
    L_b = [L_band L_tran];
    
    % Combine matrices
    %
    A_b = [A_U; -A_L];
    b = [U_b -L_b];
    
    % Set minimization vector to add up transform 
    % in transition region. For magnitude filter design
    % this will minimize the energy.  
    %
    fmin = sum(A(idx_tran,:), 1);
    
    % Call minimization routine
    %
    options = optimoptions(@linprog, 'Algorithm', 'dual-simplex', 'Display', 'off');
    
    %optimset('LargeScale', 'off', 'Algorithm', 'active-set', 'Display','off')
    
    if real_filter, 
	[x,fval,exitflag,output] = ...
            linprog(fmin, A_b, b, [],[],[],[],options);
    else
	[x,fval,exitflag,output] = ...
            linprog(fmin, A_b, b, [],[],[],[],options);
    end;

    if dbg >= 3,
	H = A * x;
	figure;
	plot_spec(f,a,d);
	[wsort, sidx] = sort(w);
	plot(w(sidx), H(sidx));
	hold on;
	plot(w(sidx), H(sidx),'rx');
	title('Frequency response calculated with A');
    end;
    
    if (exitflag == 1) % feasible
	h = fill_h(x,nhalf,real_filter, odd_filter,dbg);
	status = 'Solved'; 
    else
	h = [];
	status = 'Failed'; 
    end;
    return;

function h = fill_h(x,nhalf,real_filter,odd_filter,dbg)
% Function to fill in filter taps from optimization parameters
%
    x = x(:);
    if real_filter,
	if odd_filter, 
	    h = x(1:end);
	    h = [x(end:-1:2); h];
	else
	    h = x(1:end);
	    h = [x(end:-1:1); h];
	end;
    else
	if odd_filter, 
	    h = x(1:nhalf) + i * [0; x(nhalf+1:end)];
	    h = [conj(h(end:-1:2)); h];
	else
	    h = x(1:nhalf) + i * x(nhalf+1:end);
	    h = [conj(h(end:-1:1)); h];
	end;
    end;
    
    return;

function x0 = fill_opt_param(h0,nx,real_filter,odd_filter,lb_resp,dbg)
% Function to fill in optimization vector given a filter 
% description.  Doesn't really work well unless the previous
% filter shares the same oddness/evenness so I will ignore
% that case for now
%
    % Initialize x0 and x lengths
    %
    x0 = zeros(nx,1);
    if real_filter, 
	nx_half = nx;
    else
	if odd_filter, 
	    nx_half = (nx+1)/2;
	else
	    nx_half = nx/2;
	end;
    end;
    
    % Initialize whether FFT init should be used
    %
    fft_init = 0;
    
    if isempty(h0), 
	fft_init = 1;
    else
	% Get lengths
	%
	nh = length(h0);
	nh_half = ceil(nh/2);
    
	% Check to see if the oddness/evenness of the filter
	% is the same
	%
	if ((odd_filter && ~bitget(nh,1)) || ...
	    (~odd_filter && bitget(nh,1)))
	    fft_init = 1;
	end;
    end;
    
    if (fft_init), 
	% Get FFT of lower bound to use as initialization
	%
	nh_half = nx_half;
	if odd_filter, 
	    nh = nh_half*2 - 1;
	else
	    nh = nh_half * 2;
	end;
	h0 = fftr(hamming(length(lb_resp)).*lb_resp(:),nh);
    end;

    % Now copy taps
    %
    if odd_filter, 
	if real_filter, 
	    x0(1:min(nx_half,nh_half)) = ...
		real(h0(nh_half:nh_half+min(nx_half,nh_half)-1));
	else
	    x0(1:min(nx_half,nh_half)) = ...
		real(h0(nh_half:nh_half+min(nx_half,nh_half)-1));
	    x0(nx_half+1:nx_half+min(nx_half,nh_half)-1) = ...
		imag(h0(nh_half+1:nh_half+min(nx_half,nh_half)-1));
	end;
    else
	if real_filter, 
	    x0(1:min(nx_half,nh_half)) = ...
		real(h0(nh_half+1:nh_half+min(nx_half,nh_half)));
	else
	    x0(1:min(nx_half,nh_half)) = ...
		real(h0(nh_half+1:nh_half+min(nx_half,nh_half)));
	    x0(nx_half+1:nx_half+min(nx_half,nh_half)) = ...
		imag(h0(nh_half+1:nh_half+min(nx_half,nh_half)));
	end;
    end;
return;
    
   
    
    
    
