function [h, status] = fir_qprog_phs(n, f, ac, dc, x0, dbg)
% FIR_QPROG_PHS - FIR filter design using quadratic programming
%   
% Design n-tap linear-phase filter that meets multiband frequency specification,
% including nominal phase in each band.  
%
% function [h, status] = fir_qprog_phs(n, f, a, d, dbg)
%    
% Inputs: --- similar to cfirpm
%   n: number of taps returned
%   f: frequency bands (-1->1)
%   a: amplitude for each band  - complex  NOTE: does NOT support sloped
%       bands 
%      - magnitude describes target magnitude for band    
%      - phase describes target phase for band    
%   d: ripple  in bands - complex
%      - magnitude describes +/- magnitude ripple for band    
%      - phase describes +/- phase ripple for band    
%   dbg: flag to turn on debugging statements/plots
%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2013 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Header: /home/adam/cvsroot/src/ss/fir_qprog_phs.m,v 1.6 2013/08/15 03:18:58 adam Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    % Default value for dbg
    %
    if nargin < 6, 
	dbg = 0;
    end;
    
    % Find how many bands
    %
    nband = length(f)/2;

    % Make sure specification doesn't include sloped bands
    %
    for band = 1:nband,
        if (ac(band*2-1) ~= ac(band*2))
            error('Does not support sloped bands');
        end
    end

    % Decimate ac so it is just nband elements
    %
    a = ac(1:2:end);
    
    % Get magnitudes / phase 
    a = abs(a);
    aphs = angle(a);
    d = abs(dc);
    dphs = angle(dc);
    
    
    % Check input parameters, informing user that any stopband (e.g. a 
    % band where magnitude + ripple specs straddle 0 cannot have a 
    % phase target or ripple associated with it)
    %
    for band = 1:nband
        if ((a(band) + d(band)) * (a(band) - d(band)) < 0)
            if (a(band) ~= 0) || (dphs(band) ~= 0)
                error('Bands straddling 0 must have a = 0, angle(d) = 0');
            end
        end
    end
    
    % Check phase spec / magnitudes to see if approximation of concave
    % surface by linear segment is imposing too much error 
    %
    err_tol = 0.05;                    
    for band = 1:nband
        bandwarn = 0;
        if a(band) ~= 0
            magerr_inner = (a(band) - d(band)) * (sec(dphs(band)) - 1);
            if (magerr_inner >= 2 * d(band))
                warning('Reducing phase ripple to that feasible');
                dphs(band) = 0.99*acos((a(band) - d(band)) / (a(band)+d(band)));
            elseif ((magerr_inner > err_tol * 2 * d(band)) && (dbg > 0))
                warning(sprintf(['Band: %d, Linear approximation to magnitude spec ' ...
                                 'decreasing mag ripple by %3.1f%%'], band, 100*magerr_inner/(2*d(band))));
            end
        end
    end
    
    % Work out phases of endpoints to use for upper magnitude piecewise
    % linear segments approximation
    %
    n_phs_tran = ceil(2 * pi / acos(1-err_tol));
    amax = max(a + d);                  % Max magnitude for filter
    phs_tran = [0:n_phs_tran]/n_phs_tran * 2 * pi;
    for band = 1:nband
        if (a(band) == 0)
            phs_band{band} = [0:n_phs_tran]/n_phs_tran * 2 * pi;
        else
            phs_tol = acos(1 - (err_tol * 2 * d(band)));
            n_phs = ceil(2 * dphs(band) / phs_tol);
            phs_band{band} = ([0:n_phs]/n_phs * 2 - 1) * dphs(band) + ...
                aphs(band);

            % Add phase points from each band to transition band pwl segments
            % if band is near max magnitude
            %
            if ( (a(band) + d(band)) >= amax * (1-err_tol) )
                phs_tran = [phs_tran (aphs(band)-dphs(band)) (aphs(band)+ ...
                                                              dphs(band))];
            end
        end
    end

    % Prune and sort phase points in phs_tran
    %
    phs_tran = mod(phs_tran, 2 * pi);
    phs_tran = unique([phs_tran 0 2*pi]); % Between 0 and 2pi
    
    % Plot PWL segments for debug purposes
    %
    if (dbg > 0)
        figure;
        plot([0 (1+i)*1e-9]);
        hold on;
        
        % Plot transition limits first
        %
        tran_pts = amax * exp(i*phs_tran);
        plot(tran_pts,'k--');
        
        % Now plot bands
        %
        band_cmap = jet(nband);

        for band = 1:nband
            % Plot ideal magnitude specs
            if a(band) ~= 0,
                phs = aphs(band) + (([0:50]/50) * 2 - 1) * angle(dc(band));
            else
                phs = aphs(band) + (([0:50]/50) * 2 - 1) * pi;
            end
            plot((a(band)+d(band)) * exp(i*phs), '--',...
                 'Color', band_cmap(band,:));
            plot((a(band)-d(band)) * exp(i*phs), '--',...
                 'Color', band_cmap(band,:));
            
            % Upper magnitude is PWL approximation
            %
            plot((a(band)+d(band)) * exp(i*phs_band{band}), ...
                 'Color', band_cmap(band,:));
            
            % Lower magnitude is single line
            %
            if (a(band) ~= 0)
                pts_real = (a(band) - d(band)) * ones(1,51);
                phs = (([0:50]/50) * 2 - 1) * dphs(band);
                pts_imag = (a(band) - d(band)) * tan(phs);
                pts = (pts_real + i *pts_imag) * exp(i*aphs(band));
                plot(pts, 'Color', band_cmap(band,:));
            end
        end
        drawnow;
    end        
    
    % Scale f to -pi .. pi
    %
    f = f * pi;
    
    % Determine if filter has odd or even number of 
    % taps
    %
    if (bitget(n,1) == 1) 
	odd_filter = 1;
    else
	odd_filter = 0;
    end;
    
    % If the frequency specification has a non-zero point
    % at +/- pi, then the order must be even. A warning is 
    % printed and a failure returned if this is the case.
    %
    if (~odd_filter)
	idx = find(abs(f) == pi);
	if find(abs(ac(idx)) ~= 0)
	    warning('n odd and frequency spec non-zero at fs/2');

	    status = 'Failed';
	    h = [];
	    return;
	end;
    end;
    
    % Determine number of optimization parameters
    %
    nhalf = ceil(n/2);		% number of taps in half-side of
				% filter

    nx = 2 * n;                 % Number of optimization parameters
    
    % Create optimization arrays
    %
    oversamp = 15;
    undersamp_tran = 1;		% Undersampling factor for transition
                                % regions
    % Get first pass on w
    %
    m = 2 * oversamp * n;
    w = linspace(-pi,pi,m);
    
    % Add explicit samples to w at the edge of each specified band
    %
    w = sort([w f]);

    % Create W matrix representing DFT
    %
    if (odd_filter)
        W = exp (-i* kron(w', [-(nhalf-1):nhalf-1]));
    else
        W = exp (-i* kron(w', [-nhalf:nhalf-1]+0.5));
    end
    
    % Find indices to passbands/stopbands, and fill in upper/lower bounds
    %
    idx_band = []; 
    Au = []; Bu = [];
    Al = []; Bl = [];
    for band = 1:nband, 
	idx = find( (w >= f(band*2-1)) & (w <= f(band*2)) );
        idx_band = [idx_band idx];
        
        % Build up upper magnitude constraints
        %
        phs_diff = angle(exp(i*phs_band{band}(2)) * ...
                         exp(-i*phs_band{band}(1)));
        a_mid = (a(band) + d(band)) * cos(phs_diff/2);
        for phs_idx = 1:(length(phs_band{band})-1)
            phs_mid = phs_band{band}(phs_idx) + phs_diff/2;
            Wtmp = W(idx,:) * exp(-i*phs_mid);

            Au = [Au; real(Wtmp) -imag(Wtmp)]; % in-phase part of Wtmp * x
            Bu = [Bu; a_mid * ones(length(idx),1)];
        end
        
        % Build lower magnitude constraint for non-stopbands
        %
        if (a(band) ~= 0)
            Wtmp = W(idx,:) * exp(-i*aphs(band));
            Al = [Al; real(Wtmp) -imag(Wtmp)]; % in-phase part of Wtmp * x
            Bl = [Bl; (a(band)-d(band)) * ones(length(idx),1)];

            % Build upper phase constraint
            %
            Wtmp = W(idx,:) * exp(-i*phs_band{band}(end));
            Au = [Au; imag(Wtmp) real(Wtmp)]; % quadrature part of Wtmp * x
            Bu = [Bu; zeros(length(idx),1)];
        
            % Build lower phase constraint
            %
            Wtmp = W(idx,:) * exp(-i*phs_band{band}(1));
            Al = [Al; imag(Wtmp) real(Wtmp)]; % quadrature part of Wtmp * x
            Bl = [Bl; zeros(length(idx),1)];
        end
    end

    % Get transition indices
    %
    idx_tmp = ones(1,length(w));
    idx_tmp(idx_band) = 0;
    idx_tran = find(idx_tmp == 1);

    % Decimate w in transition regions
    %
    idx_tran = idx_tran(1:undersamp_tran:end);

    % Find weighting
    %
    wband_mtx = w(idx_band)' * ones(1,length(idx_tran));
    w_mtx = ones(length(idx_band),1) * w(idx_tran);
    dtmp = abs(angle(exp(i*(wband_mtx - w_mtx))));
    dtmp = min(dtmp);                   
    wt_tran_energy = min(4*pi/n, dtmp) / (4*pi/n);
    
    if dbg >= 3, 
	figure;
        hold on;
	plot(w(idx_band),ones(length(idx_band),1),'*');
        plot(w(idx_tran), wt_tran_energy, 'g+');
        drawnow;
    end;
    
    % Build up transition band magnitude constraints
    %
    limit_tran = 1;
    if limit_tran
        for idx = 1:(length(phs_tran)-1)
            phs_diff = phs_tran(idx+1)-phs_tran(idx);
            phs_mid = phs_tran(idx)+phs_diff/2;
            Wtmp = W(idx_tran,:) * exp(-i*phs_mid);

            Au = [Au; real(Wtmp) -imag(Wtmp)]; % in-phase part of Wtmp * x
            Bu = [Bu; amax * cos(phs_diff/2) * ones(length(idx_tran),1)];
        end
    end
    
    % Combine matrices
    %
    A = [Au; -Al];
    B = [Bu; -Bl];
    
    % Build matrix H to minimize energy of filter
    %
    minimize_total_energy = 1;
    if minimize_total_energy
        H = eye(nx);                        % Minimize total energy
        H = sparse(H);
    else
        Ar_tran = [real(W(idx_tran,:)) -imag(W(idx_tran,:))];
        Ai_tran = [imag(W(idx_tran,:)) real(W(idx_tran,:))];
        Ari_tran = diag([sqrt(wt_tran_energy) sqrt(wt_tran_energy)]) * [Ar_tran; Ai_tran];
        H = Ari_tran' * Ari_tran;
    end

    fmin = zeros(1,nx);
    
    % Call minimization routine
    %
    x0 = [];
    [x,fval,exitflag,output] = ...
        quadprog(H, fmin, A, B, [],[],[],[],x0,...
                            optimset('Algorithm', 'interior-point-convex', ...
                                     'Display','off'));
    %                     optimset('Algorithm', 'active-set', 'Display','iter', ...
    %                          'MaxIter', 500));

    
    % exitflag

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
        h = x(1:n) + i * x(n+1:end);
	H = W * h;
	[wsort, sidx] = sort(w);
        hsort = H(sidx);
	figure;
	plot_spec_phs(f,ac,dc,wsort,hsort);
	title('Frequency response calculated with W');
        figure;
        plot(abs(h));
        drawnow;
    end;
    
    if (exitflag == 1) % feasible
	h = x(1:n) + i * x(n+1:end);
	status = 'Solved'; 
    else
	h = [];
	status = 'Failed'; 
    end;
    
    return;

    
    
    
