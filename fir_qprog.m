
function [h, status] = fir_qprog(n, f, a, d, dbg)
% FIR_QPROG - FIR filter design using quadratic programming
%   
% Design n-tap linear-phase filter that meets multiband frequency 
% specification.  
%
% Also constrain H(f) > max(0, min(a(1:2:end)-d), min(a(2:2:end)-d))
% Also minimizes sum(|H(f)|^2) in transition bands.
%
% function [h, status] = fir_qprog(n, f, a, d, dbg)
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
% $Header: /home/adam/cvsroot/src/ss/fir_qprog.m,v 1.2 2013/08/15 17:10:31 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    % Default value for dbg
    %
    if nargin < 5, 
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
    
    % Set H to minimize total energy in filter
    % Set fmin to 0
    H = eye(nx);
    fmin = zeros(1,nx);
    
    % Call minimization routine
    %
    x0 = [];
    if real_filter, 
	[x,fval,exitflag,output] = ...
            quadprog(H, fmin, A_b, b, [],[],[],[],x0,...
                    optimset('Algorithm', 'interior-point-convex', ...
                             'Display','off'));
    else
	[x,fval,exitflag,output] = ...
            quadprog(H, fmin, A_b, b, [],[],[],[],x0,...
                    optimset('LargeScale','off', 'Algorithm', 'interior-point-convex', 'Display','off'));
    end;

    if dbg >= 2,
        fprintf(1,'Exitflag: %d\n', exitflag);
        switch(exitflag)
          case 1
            fprintf(1,'First order optimality conditions satisfied\n');
          case 0 
            fprintf(1,'Maximum number of iterations exceeded\n');
          case -2 
            fprintf(1,'No feasible point found\n');
          case -3 
            fprintf(1,'Problem is unbounded\n');
          case -6 
            fprintf(1,'Non-convex problem detected\n');
          case 3 
            fprintf(1,'Change in objective function too small\n');
          case -4 
            fprintf(1,['Current search direction is not a descent direction; ' ...
                       'no further progress can be made.\n']);
          case 4 
            fprintf(1,'Local minimizer found\n');
          case -7 
            fprintf(1,['Magnitude of search direction became too small; no ' ...
                       'further progress can be made. The problem is ill-posed ' ...
                       'or badly conditioned.\n']);
          otherwise 
            fprintf(1,'Exitflag not recognized\n');
        end

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

