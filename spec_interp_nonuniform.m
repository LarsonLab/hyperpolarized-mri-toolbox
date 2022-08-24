function [hi, t_all] = spec_interp_nonuniform(h, ni, off, f, dbg)
% SPEC_INTERP_NONUNIFORM - Interpolate filter to nonuniform taps while
% keeping spectral response consistent
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


if isempty(off)
    off = 0.5;
end

  N = length(h);
    mult_factor = 15;
    w = linspace(-pi, pi, 2*mult_factor*N);
    
    % Try weighting samples that are in frequency band 
    % more heavily
    %
    f = f * pi;
    idx_band = [];
    nband = length(f)/2;
    for band = 1:nband, 
	idx = find( (w >= f(band*2-1)) & (w <= f(band*2)) );
	idx_band = [idx_band idx];
    end;
    wt_band = 10;
    wt = ones(length(w),1);
    wt(idx_band) = wt_band;
    
    % Get reference transform
    %
    t_ref = [0:N-1];
    Wref = exp(-i*kron(w', t_ref));
    Fref = Wref * h(:);
    Fref_wt = wt .* Fref;
    
    if (dbg >= 2), 
	filt_fig = figure;
    end;
 
    t_off = ([1:ni]-(ni+1)/2)/(ni-1)*2*off;
    t_all = zeros(ni, N);
    
    for idx = [1:ni]
	% Get actual sampling positions
	%
        t_all(idx,:) = t_ref + t_off(idx)* (-1).^[0:N-1];
        t_act = t_all(idx,:);
	
	if (dbg >= 2)
	    figure(filt_fig);
	    subplot(411);
	    stem([t_ref.' t_act.'], ones(length(t_ref),2));
	    legend('Reference', 'Actual');
	    title('Sampling Locations');
	end;

	% Get actual, add weights
	%
	Wact = exp(-i*kron(w', t_act));
	Wact_wt = repmat(wt,1,length(t_act)) .* Wact;
        
	% Get new filter
	%
	hi(idx,:) = pinv(Wact_wt)*Fref_wt;
    
    	if ( (dbg >= 2) && (rem(idx,5)==0) ), 
	    plot_db = 0;
            figure(filt_fig);
            
            Fact = Wact * h(:);
            Fact_fix = Wact * hi(idx,:).';

            subplot(4,1,2);
            hold off;
            plot(abs(h));
            hold on;
            plot(abs(hi(idx)), 'r--');
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

	    fprintf(1,'Offset: %f   -- Hit any key to continue\n', t_off(idx));
	    pause;
	end;

    
    end;
    
    hi = hi(:);
    t_all = t_all(:);
    
    return;
    
    
    


