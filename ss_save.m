function ss_save(g,rf,ang,thk, isodelay, format, fspec, a_angs)
% SS_SAVE - Save spectral-spatial pulse 
% Uses Chuck Cunningham's format for GE systems, and creates associated
% .dat-file
% Pulse parameters saved in header for Varian fules
%   
%  ss_save(g,rf,ang,thk, isodelay, format, fspec, a_angs)
%
%  g - in G/cm    
%  rf - in G
%  ang - flip angle in radians
%  thk - thickness in cm
%  isodelay (optional) - delay from in-phase point to end of pulse (GE
%  definition)
%  format (optional) - 'GE' (default), 'Varian'  
%  fspec (optional) - frequency bands (Hz) to write in file
%  a_angs (optional) - band amplidutes (radians) to write in file

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


ss_globals;
    
if (nargin < 5) || isempty(isodelay)
    isodelay = length(rf)*SS_TS/ 2;
end

if (nargin < 6) || isempty(format)
    format = 'GE';
end

    switch format,
        case 'GE'
            % force even number of samples due to some GE sequences having
            % issues loading odd number of samples
            if rem(length(rf))
                rf = [rf(:), 0];
                g = [g(:), 0];
            end
            
        case 'Varian'
            
        otherwise
            error(sprintf(['Format save type of: %s not' ...
                ' recognized'], format));
    end;
    
    maxg = max(abs(g));
    if maxg ~= 0, 
	gn = g/maxg;
    else
	gn = g;
    end;
    
    maxrf = max(abs(rf));
    rfn = rf / maxrf;

    root_fname = input('Root file name: (leave empty to not save) ', 's');
    if isempty(root_fname)
	fprintf(1,'Not saving files \n');
	return;
    end;

    % calculate pulse parameters
    nrf = length(rf);
    abswidth = sum(abs(rfn))/nrf;
    effwidth = sum(abs(rfn).^2)/nrf;
    area = sum(abs(rfn))/nrf;
    pon = (rfn >= 0.00001);
    temp_pw = 0;
    max_pw = 0;
    for n=1:nrf
	temp_pw = temp_pw + pon(n);
	if (and(pon(n) == 0, temp_pw ~= 0))
	    max_pw = max(max_pw, temp_pw);
	    temp_pw = 0;
	end;
    end;
    max_pw = max_pw / n;
    
    dty_cyc = sum(abs(rfn) > 0.2236)/nrf; 
    if dty_cyc < max_pw, 
	dty_cyc = max_pw;
    end;

    max_b1 = max(abs(rf));
    int_b1_sqr = sum(abs(rf).^2 * SS_TS * 1e3);
    rms_b1 = sqrt(sum(abs(rf).^2))/nrf;
    thk_scale = thk * SS_GAMMA / SS_GAMMA_HYDROGEN * 10;
    
    
    % Allow magnitude of RF to go negative, this will
    % help reduce sensitivity to theta modulation since
    % there is a possible delay on the system for this case
    %
    
    % Find pi jumps in phase and remove
    % Do this by doubling phase and unwrapping 2*pi jumps
    %
    dang_rf = 2*angle(rfn);
    dang_rf = dang_rf - dang_rf(1);
    dang_rf = unwrap(dang_rf, 0.98*pi);
    ang_rf = mod(dang_rf/2+pi,2*pi)-pi;
    
    mag_rf = real(rfn .* exp(-i*ang_rf));
    
    if 1, 
	figure;
	subplot(211)
	t = 1:length(rf);
	plot(abs(rf));
	hold on;
	plot(maxrf*mag_rf,'r--')
	
	subplot(212);
	plot(real(rf));
	hold on;
	plot(imag(rf), 'b--');
	plot(real(maxrf*mag_rf.*exp(i*ang_rf)), 'r:');
	plot(imag(maxrf*mag_rf.*exp(i*ang_rf)), 'r:');
    end;

    
    
    switch (format)
        case 'GE'
    
            dat_name = sprintf('%s.dat', root_fname);
            fid = fopen(dat_name, 'w');
            if fid == -1, 
            fprintf(1, 'Error opening %s \n', dat_name);
            return;
            end;

            fprintf(fid,'%10d \t\t #Spectral-Spatial\n', 1);
            fprintf(fid,'%10d \t\t #res\n', length(rf));
            fprintf(fid,'%10d \t\t #pw\n',round(length(rf)*SS_TS*1e6));
            fprintf(fid,'%10.7f \t\t #nom_flip \n',ang*180/pi);
            fprintf(fid,'%10.7f \t\t #abswidth \n',abswidth);
            fprintf(fid,'%10.7f \t\t #effwidth \n',effwidth);
            fprintf(fid,'%10.7f \t\t #area \n',area);
            fprintf(fid,'%10.7f \t\t #dtycyc \n',dty_cyc);
            fprintf(fid,'%10.7f \t\t #maxpw \n',max_pw);
            gamscale = SS_GAMMA/SS_GAMMA_HYDROGEN; % GE assumes max B1 is for application at 1H
            fprintf(fid,'%10.7f \t\t #max_b1 \n',max_b1 * gamscale); 
            fprintf(fid,'%10.7f \t\t #max_int_b1_sqr \n',int_b1_sqr* gamscale^2);
            fprintf(fid,'%10.7f \t\t #max_rms_b1 \n',rms_b1* gamscale^2);
            fprintf(fid,'%10.3f \t\t #a_gzs \n',maxg);
            fprintf(fid,'%10.3f \t\t #nom_thk(mm) \n',thk * gamscale * 10);
            fprintf(fid,'%10d \t\t #isodelay\n',round(isodelay*1e6));
            fprintf(fid,'%10d \t\t #g_pow \n',0);
            fprintf(fid,'%10d \t\t #g_pos_pow \n',0);
            fprintf(fid,'%10d \t\t #g_neg_pow \n',0);
            fprintf(fid,'%10d \t\t #g_abs \n',0);
            fprintf(fid,'%10d \t\t #g_dgdt \n',0);
            fprintf(fid,'%10d \t\t #g_pwm \n',0);
            fprintf(fid,'%10d \t\t #g_pwm_abs \n',0);
            fprintf(fid,'# *************************************\n');
            if (nargin == 7)
                for b = 1:length(a_angs)
                   fprintf(fid,'# Band %d: [%.2f, %.2f] Hz, %.2f degree flip\n', ...
                       b, fspec(2*b-1), fspec(2*b), a_angs(b)*180/pi);
                end
            end
            
            fclose(fid);

            % Now write out RF and Gradient 
            %
            rho_fname = sprintf('%s.rho', root_fname);
            signa(mag_rf,rho_fname,1);

            theta_fname = sprintf('%s.pha', root_fname);
            signa(-ang_rf,theta_fname,1/pi);

            g_fname = sprintf('%s.grd', root_fname);
            signa(gn,g_fname,1);
    
        case 'Varian'
            
            fid = fopen(sprintf('%s.RF', root_fname),'wt');
            fprintf(fid,'# %s\n', sprintf('%s.RF', root_fname));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# Spectral-spatial Matlab Package\n');
            fprintf(fid,'# Nucleus = %s\n',SS_NUCLEUS);
            fprintf(fid,'# Duration = %d us\n',round(length(rf)*SS_TS*1e6));
            fprintf(fid,'# Isodelay = %d us\n',round(isodelay*1e6));
            fprintf(fid,'# Resolution = %d us\n', SS_TS*1e6);
            fprintf(fid,'# Flip = %.2f degrees\n',ang*180/pi);
            fprintf(fid,'# Max B1 = %.4f Gauss\n',max_b1);
             if (nargin == 7)
                for b = 1:length(a_angs)
                   fprintf(fid,'# Band %d: [%.2f, %.2f] Hz, %.2f degree flip\n', ...
                       b, fspec(2*b-1), fspec(2*b), a_angs(b)*180/pi);
                end
            end
           fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# VERSION       Matlab\n');
            fprintf(fid,'# TYPE          selective\n');
            fprintf(fid,'# MODULATION    amplitude\n');
            fprintf(fid,'# EXCITEWIDTH   -1.0000\n');
            fprintf(fid,'# INVERTWIDTH   -1.0000\n');  
            fprintf(fid,'# INTEGRAL      %1.4f\n',    sum(abs(rfn))/length(rfn)); 
            fprintf(fid,'# RF_FRACTION   -1.0000\n');
            fprintf(fid,'# STEPS         %d\n',length(rfn));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'%3.2f  %4.2f  1.0\n',[-angle(rfn(:)).'*180/pi; 1023*abs(rfn(:)).']);    
            fclose(fid);

            fid = fopen(sprintf('%s.GRD', root_fname),'wt');
            fprintf(fid,'# %s\n', sprintf('%s.GRD', root_fname));
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'# Spectral-spatial Matlab Package\n');
            fprintf(fid,'# Nucleus = %s\n',SS_NUCLEUS);
            fprintf(fid,'# Duration = %d us\n',round(length(g)*SS_TS*1e6));
            fprintf(fid,'# Resolution = %d us\n', SS_TS*1e6);
            fprintf(fid,'# Points = %d\n',length(g));
            fprintf(fid,'# Max Gradient Strength = %.4f Gauss/cm\n',maxg);
            fprintf(fid,'# (Max Gradient Strength Constraint = %.2f Gauss/cm)\n',SS_MXG);
            fprintf(fid,'# Max Slew Rate = %.2f Gauss/cm/ms\n',SS_MXS);
            fprintf(fid,'# Slice thickness = %.1f mm\n',thk*10);
            fprintf(fid,'# ***************************************************\n');
            fprintf(fid,'%d  1\n',round(32767*gn));
            fclose(fid);
    
    end
    
    
    
    

    
