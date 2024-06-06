function [f_a, a_a, d_a, f_off] = ss_alias(f, a, d, f_off, fs, sym, fig, threshold_edge)
% SS_ALIAS - Alias filter specifications into effective bandwidth
%   
% function [f_a, a_a, d_a, f_off] = ss_alias(f, a, d, f_off, fs, sym, fig, threshold_edge)
%
% Inputs:     
% f - frequency band edges, in Hz, each band monotonically increasing, 
%     length(f) = 2 * number of bands
% a - amplitude for each band [0..1]
% d - ripple weighting for each band 
% f_off - offset frequency in Hz (if [], then will be 
%         set to include as many passbands within one sampling
%         interval as possible    
% fs - sampling frequency (bandwidth), in Hz
% sym - flag indicating whether frequency response should be symmetric  
% threshold_edge - [optional], threshold of the min distance to the 
%                  normalized spectrum edge [-1,1], usually 0.15-0.2
% fig - [optional], figure to plot
%
% Outputs: 
% f_a - Aliased frequency band edges into normalized freq [-1..1]
% a_a - amplitude of aliased bands
% d_a - ripple weighting of aliased bands
% f_off - offset frequency in Hz
%     
% Modification by Hong Shang, June 2013
% 1 Display the aliased filter spec if with one more input, fig, plot 
% normalized band before and after aliasing and compatible flag, to help
% visualize the choice of BW 
% 2 the choice of f_off when sym=0, for both including as many bands in the
% main spectrum and getting rid of bands at spectrum edge.
% 3 ripple can be different if overlap
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
% $Header: /home/adam/cvsroot/src/ss/ss_alias.m,v 1.12 2014/05/22 20:43:59 peder Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check input parameters
    %
    if nargin < 6, 
        error('Usage: [f_a, a_a, d_a, f_off] = ss_alias(f, a, d, f_off, fs, sym)\n');
    elseif nargin == 6
        threshold_edge = 0.03;
        display_flag = 0;
    elseif nargin == 7
        if ishandle(fig)
            display_flag = 1;
            threshold_edge = 0.03;
        else
            display_flag = 0;
            threshold_edge = fig;
        end;
    elseif nargin == 8
        display_flag = 1;
    end;
        
    if mod(length(f), 2) ~= 0, 
        error('Frequency vector (f) of band pairs must be even length');
    end;
    
    nf = length(f)/2;
    if length(d) ~= nf, 
        error('Ripple vector incorrect length');
    end;
    if length(a) ~= nf, 
        error('Amplitude vector incorrect length');
    end;
    
    df = diff(f);
    if (min(df(1:2:end)) <= 0) 
        error('Frequency bands must be monotonic increasing');
    end;

    if (max(df(1:2:end)) >= fs)  % there is one band longer than bandwidth, impossible to design such filter
        f_a = [];  a_a = [];  d_a = [];  f_off = 0;
        return;
    end;
    
    % If sym flag then either set or check f_off
    if (sym), 
        if ~isempty(f_off), 
	    % Check specified frequency offset to make sure it is < midpoint 
	    % of first band or > than midpoint of top band
            if ((f_off > (f(2)-f(1))/2) &&  (f_off < (f(end) - f(end-1))/2) )
                error('f_off not at band edges');
            end;
        else
	    % Set f_off to left band edge
            f_off = (f(2)+f(1))/2;
        end;
    else
        % Set f_off to include as many bands as possible in one sampling interval
        if isempty(f_off)
            if (max(f) - min(f)) <= fs,   % all band can exist in the center spectrum
                f_off = (max(f) + min(f))/2;
            else
            % Find how many bands are in -fs/2..fs/2
%                 f_l = f(1:2:end);
%                 f_u = f(2:2:end);
%                 f_off_test = linspace(min(f)+fs/2, max(f)-fs/2, 100);
%                 band_in = zeros(length(f_off_test),length(f_l));
%                 nband = zeros(1,length(f_off_test));
%                 for off_idx = 1:length(f_off_test), 
%                     f_off = f_off_test(off_idx);
%                     band_in(off_idx,:) = ((f_l >= f_off-fs/2) & (f_l <= f_off+fs/2) &  (f_u >= f_off-fs/2) & (f_u <= f_off+fs/2));
%                     nband(off_idx) = sum(band_in(off_idx,:));
%                 end;
%                 [~, off_idx] = max(nband);
%                 in_idx = find(band_in(off_idx,:));
%                 f_off = (f(in_idx(1)*2-1) + f(in_idx(end)*2))/2;
                
                % modification of choosing f_off
                f_l = f(1:2:end);
                f_u = f(2:2:end);
                f_off_test = linspace(min(f), max(f), 500);  % Hz
                band_in = zeros(length(f_off_test),length(f_l));
                nband = zeros(1,length(f_off_test));
                edge_distance = zeros(1,length(f_off_test));
                split_flag = zeros(1,length(f_off_test));
                for off_idx = 1:length(f_off_test), 
                    f_off = f_off_test(off_idx);
                    band_in(off_idx,:) = ((f_l >= f_off-fs/2) & (f_l <= f_off+fs/2) &  (f_u >= f_off-fs/2) & (f_u <= f_off+fs/2));
                    nband(off_idx) = sum(band_in(off_idx,:));
                    
                    % pre-calculation of Normalized and aliased frequency
                    fnorm_pre = zeros(size(f));
                    for idx = 1:nf
                        % Normalized frequencies
                        fa1 = (f(idx*2-1) - f_off) / (fs/2); 
                        fa2 = (f(idx*2) - f_off) / (fs/2); 
                        % Get aliased frequencies
                        fa1 = mod(fa1+1,2)-1;
                        fa2 = mod(fa2+1,2)-1;
                        % get rid of confusing about whether at -1, or 1
                        if (fa2 < fa1) && (fa2 == -1)
                            fa2 = 1;
                        end;
                        if (fa2 < fa1) && (fa1 == 1)
                            fa1 = -1;
                        end;
                        % check whether a band is split
                        if  fa2 <  fa1  
                            split_flag(off_idx) = 1;
                        end
                        fnorm_pre(idx*2-1) = fa1;
                        fnorm_pre(idx*2) = fa2;                
                    end
                    edge_distance(off_idx) = min( [min(fnorm_pre(1:2:end)+1), min(1-fnorm_pre(2:2:end))] );
                    if split_flag(off_idx) > 0.5;  edge_distance(off_idx) = 0; end;
                end;
                % figure; 
                % subplot(3,1,1); plot(f_off_test,edge_distance,'b.'); xlabel('offset frequency, Hz'); title('min distance to the edge'); axis tight; set(gca,'ylim',[-0.1 0.7]);
                % subplot(3,1,2); plot(f_off_test,split_flag,'b.');    xlabel('offset frequency, Hz'); title('whether some bands are split');  axis tight; set(gca,'ylim',[0 1.5]);
                % subplot(3,1,3); plot(f_off_test,nband,'b.');         xlabel('offset frequency, Hz'); title('how many bands in the main spectrum');  axis tight; set(gca,'ylim',[0 nf]);
                
                % find the optimal f_off
                idx_1 = find( split_flag < 0.5 );
                %f_off_test = f_off_test(idx_1); 
                %edge_distance = edge_distance(idx_1); 
                %nband = nband(idx_1);
                idx_2 = find( nband == max(nband) );               
                [value_3,idx_3] = max(edge_distance(idx_2));       
                if value_3 > threshold_edge
                    f_off = f_off_test(idx_2(idx_3));
                else
                    idx_4 = find( nband >= (max(nband)-1) );
                    [~,idx_5] = max(edge_distance(idx_4));
                    f_off = f_off_test(idx_4(idx_5));
                end;
            end;
        end;
    end;
	
    % If symmetric frequency response, then construct mirrored response, update number of bands
    if (sym)
        f_h = f - f_off;
        f = [-f_h(end:-1:1) f_h] + f_off;
        a = [a(end:-1:1) a];
        d = [d(end:-1:1) d];
        nf = length(f)/2;
    end;

    % Go through each band pair, determine aliased frequencies 
    % in normalized space, and split if necessary
    new_idx = 1;
    fnorm = zeros(size(f));
    for idx = 1:nf, 
        % Normalize frequencies
        fnorm(idx*2-1) = (f(idx*2-1) - f_off) / (fs/2); 
        fnorm(idx*2) = (f(idx*2) - f_off) / (fs/2); 

        % Get aliased frequencies
        fa1 = mod(fnorm(idx*2-1)+1,2)-1;
        fa2 = mod(fnorm(idx*2)+1,2)-1;
	
        % in case one band is interrupted into two ends of main spectrum
        % Check to see if endpoints can be shifted to form one single band
        if (fa2 < fa1) && (fa2 == -1)
            fa2 = 1;
        end;
        if (fa2 < fa1) && (fa1 == 1)
            fa1 = -1;
        end;
	
        if (fa2 < fa1),  % if still split, then create a new band
            f_a(new_idx*2-1) = -1;
            f_a(new_idx*2) = fa2;
            a_a(new_idx) = a(idx);
            d_a(new_idx) = d(idx);
            new_idx = new_idx + 1;

            f_a(new_idx*2-1) = fa1;
            f_a(new_idx*2) = 1;
            a_a(new_idx) = a(idx);
            d_a(new_idx) = d(idx);
            new_idx = new_idx + 1;
        else
            f_a(new_idx*2-1) = fa1;
            f_a(new_idx*2) = fa2;
            a_a(new_idx) = a(idx);
            d_a(new_idx) = d(idx);
            new_idx = new_idx + 1;
        end;
    end;

    % additional display 
    if display_flag == 1
        clf(fig);  figure(fig); 
        subplot(2,1,1); hold on;
        for i = 1:length(a)
            f_sub = [fnorm(2*i-1),fnorm(2*i)];
            bound_down = (a(i) - d(i))*ones(size(f_sub)) ;
            bound_up   = (a(i) + d(i))*ones(size(f_sub));
            plot(f_sub,bound_down,'r-',f_sub,bound_up,'r-');
        end; 
        hold off;  xlabel('normalized frequency'); title(['filter spec before aliasing at fs = ',num2str(fs),'Hz']); set(gca,'xlim',[min(fnorm)-0.1,max(fnorm)+0.1]); set(gca,'ylim',[0,max(a)+max(d)+3e-2]);
        
        figure(fig); subplot(2,1,2); hold on;
        for i = 1:length(a_a)
            f_sub = [f_a(2*i-1),f_a(2*i)];
            bound_down = (a_a(i) - d_a(i))*ones(size(f_sub)) ;
            bound_up   = (a_a(i) + d_a(i))*ones(size(f_sub));
            plot(f_sub,bound_down,'r-',f_sub,bound_up,'r-');
        end; 
        hold off;  xlabel('normalized frequency');  title(['filter spec after aliasing at fs = ',num2str(fs),'Hz']);  set(gca,'xlim',[-1,1]); set(gca,'ylim',[0,max(a_a)+max(d_a)+3e-2]);
    end
     
    % Check for overlaps
    %   -- a little complicated
    overlap_found = 1;
    incompatible_overlap = 0;
    while (overlap_found && ~incompatible_overlap && length(f_a)/2 > 1) 
	% Check each frequency band for overlap
	nf = length(f_a)/2;
	f_tmp = [];
	a_tmp = [];
	d_tmp = [];
	for chk_idx = 1:nf, 
	    % Copy frequency, ripple, amplitude
	    f_tmp(chk_idx*2-1) = f_a(chk_idx*2-1);
	    f_tmp(chk_idx*2) = f_a(chk_idx*2);
	    d_tmp(chk_idx) = d_a(chk_idx);
	    a_tmp(chk_idx) = a_a(chk_idx);

	    % Check current band edges with all remaining bands for overlap
	    fl_1 = f_a(chk_idx*2-1);
	    fu_1 = f_a(chk_idx*2);
	    for rem_idx = chk_idx+1:nf, 
            % Initialize flags
            overlap_ok = 0;
            overlap_found = 0;
            
            % Check to see if overlap would be compatible i.e. same ripple, same amplitude
            % if ((d_a(chk_idx) == d_a(rem_idx)) && (a_a(chk_idx) == a_a(rem_idx))); overlap_ok = 1;  end;
            % only amplitude need to be consistent, ripple can be
            % different, use the smaller ripple if overlap
            if ( a_a(chk_idx) == a_a(rem_idx) ); overlap_ok = 1;  end; 
            
            % Check if band edges overlap
            fl_2 = f_a(rem_idx*2-1);
            fu_2 = f_a(rem_idx*2);
            if (fl_1 <= fl_2) && (fu_1 >= fl_2)
                % Fix overlapping band
                f_tmp(chk_idx*2-1) = fl_1;
                f_tmp(chk_idx*2) = max(fu_1, fu_2);
                d_tmp(chk_idx) = min([d_a(chk_idx),d_a(rem_idx)]);
                
                % Copy other bands
                for cp_idx = (chk_idx+1):nf,
                    if cp_idx ~= rem_idx,
                        f_tmp = [f_tmp f_a(cp_idx*2-1) f_a(cp_idx*2)];
                        a_tmp = [a_tmp a_a(cp_idx)];
                        d_tmp = [d_tmp d_a(cp_idx)];
                    end;
                end;
                overlap_found = 1;
                break;		% Start over
            elseif (fu_1 >= fu_2) && (fl_1 <= fu_2)
                % Fix overlapping bands
                f_tmp(chk_idx*2-1) = min(fl_1, fl_2);
                f_tmp(chk_idx*2) = fu_1;
                d_tmp(chk_idx) = min([d_a(chk_idx),d_a(rem_idx)]);
                
                % Copy other bands
                for cp_idx = (chk_idx+1):nf,
                    if cp_idx ~= rem_idx,
                        f_tmp = [f_tmp f_a(cp_idx*2-1) f_a(cp_idx*2)];
                        a_tmp = [a_tmp a_a(cp_idx)];
                        d_tmp = [d_tmp d_a(cp_idx)];
                    end;
                end;
                overlap_found = 1;
                break;		% Start over
            elseif (fl_1 >= fl_2) && (fu_1 <= fu_2)
                % Fix overlapping bands
                f_tmp(chk_idx*2-1) = fl_2;
                f_tmp(chk_idx*2) = fu_2;
                d_tmp(chk_idx) = min([d_a(chk_idx),d_a(rem_idx)]);
                
                % Copy other bands
                for cp_idx = (chk_idx+1):nf,
                    if cp_idx ~= rem_idx,
                        f_tmp = [f_tmp f_a(cp_idx*2-1) f_a(cp_idx*2)];
                        a_tmp = [a_tmp a_a(cp_idx)];
                        d_tmp = [d_tmp d_a(cp_idx)];
                    end;
                end;
                overlap_found = 1;
                break;		% Start over
            end;
	    end;
        if overlap_found,
            incompatible_overlap = ~overlap_ok;
            f_a = f_tmp;
            d_a = d_tmp;
            a_a = a_tmp;
            break;			% Start over from beginning
        end;
	end;
    end;
     
    % If incompatible overlap, then return empty matrices
    if incompatible_overlap, 
        f_a = []; d_a = []; a_a = [];
        if display_flag == 1    %modification with additional display
            figure(fig); subplot(2,1,2); legend('incompatible overlap');
        end
        return;
    end;
    
    % Sort frequencies into ascending order
    [f_a, idx_sort] = sort(f_a);
    a_a = a_a(idx_sort(2:2:end)/2);
    d_a = d_a(idx_sort(2:2:end)/2);

    % If symmetric response then return only positive bands
    if (sym)
	% Find right band edges that are positive, and fix return variables	
        idx_pos = find(f_a(2:2:end) > 0);
        f_tmp = [];
        for idx = 1:length(idx_pos), 
            f_tmp(idx*2-1) = f_a(idx_pos(idx)*2-1);
            f_tmp(idx*2) = f_a(idx_pos(idx)*2);
        end;
        f_a = max(0,f_tmp);
        a_a = a_a(idx_pos);
        d_a = d_a(idx_pos);
    end;
    
    return;	
	
    
