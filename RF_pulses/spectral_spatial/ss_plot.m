function [f,z,m] = ss_plot(g, rf, samp, ptype, fov, bw, gamma, fplot, isodelay)
% SS_PLOT - Plot performance of spectral-spatial
%   
%  function [f,z,m] = ss_plot(g, rf, samp, [ptype], [fov], [bw], [gamma],
%  [fplot], [isodelay])
%
% INPUTS
%    g - gradient in G/cm 
%    rf - RF in G    
%    samp - sample period in s
%    [ptype] - pulse type: 'ex', 'se', 'sat', 'inv'
%    [fov] - Spatial fov in cm to plot
%    [bw] - Spectral bw in Hz to plot
%    [gamma] - Gamma to be used, Default:4257
%    [fplot] - frequencies to plot spatial profiles
%    [isodelay] - Unwind spectral phase shift for given isodelay - default: 0
%
% OUTPUTS
%    f - frequency points plotted (Hz)
%    z - spatial points plotted (cm)
%    m - simulated magnetization
%
%    If optional parameters [] are set to empty matrices, 
%    it attempts to estimate reasonable values
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


    % Some defines
    %
    NPLOT = 100;			% Number plot points
    
    if nargin < 3, 
	error(['Usage: ss_plot(g, rf, samp, [ptype], [fov], [bw], [gamma],' ...
	       ' [fp], [fs])']);
    end;
    
    if ~isreal(g), 
	error('Only works with 1-D gradient');
    end;
    
    if length(g) ~= length(rf), 
	error('RF and gradient must be same length');
    end;
    g = g(:);
    rf = rf(:);

    % Check ptype
    if ( (nargin < 4) || isempty(ptype) ), 
    ptype = 'ex';
    else
      switch ptype, 
        case {'ex', 'se', 'sat', 'inv'}
        otherwise
          error(sprintf(['SS_PLOT: pulse type (ptype) of: %s not' ...
			 ' recognized'], ptype));
      end;
    end
    
    if ( (nargin < 5) || isempty(fov)) 
        fov = 4.0;
    end
    
    if ( (nargin < 6) || isempty(bw) ), 
	% Take Fourier transform of gradient and find peak 
	% component
	g_fft = fft(g);
	res = length(g);
	res_d2 = fix(res/2);
	[mxval, idx] = max(abs(g_fft(1:res_d2)));
	freq_peak = 1/samp * idx/length(g);
	
	% Plot will be from -bw/2:bw/2
	%
	bwidth = round(4 * freq_peak / 100) * 100;
	bw = [-bwidth/2 bwidth/2];
    elseif length(bw) == 1
        bw = [-bw/2 bw/2];
    end
    
    if ( (nargin < 5) || isempty(fov) ), 
	% Estimate FOV to plot from k-space swing
	%
	kz = cumsum(gamma * g * samp);
	kz_peak = (max(kz) - min(kz))/2;
	fov = ceil(16 / kz_peak);
    end;
    
    if ( (nargin < 7) || isempty(gamma) ), 
	gamma = 4257;
    end;

    if ( (nargin < 9) || isempty(isodelay) ), 
	isodelay = 0;
    end;
    
    
    

    % Get x and y vectors for abr
    %
    dbw = diff(bw) / (NPLOT-1);
    f = [bw(1):dbw:bw(2)];

    dfov = fov / (NPLOT-1);
    z = [-fov/2:dfov:fov/2];

    % Convert RF to a rotation in radians
    %
    rf_rot = 2 * pi * gamma * rf(:) * samp;

    % Build gradient that gives rotation  
    % in radians when scaled by "f" (off-resonance in Hz) 
    %
    gf_rot = 2 * pi * samp * ones(size(g(:)));
    
    % Convert gradient to a rotation in radians
    % when scaled by "z"
    %
    gz_rot = 2 * pi * gamma * g(:) * samp;
    
    % Add single samples at end of all waveforms to account for isodelay correction
    %
    rf_rot = [rf_rot;0];
    gf_rot = [gf_rot;(-2 * pi * isodelay)];
    gz_rot = [gz_rot;0];
    
    % Get Mxy now
    %
    m = calc_mag(rf_rot, gz_rot+i*gf_rot, z, f, ptype);
    
    % Make plots now 
    %
    figure; 
    
    % RF
    %
    t = [0:length(g)-1] * samp * 1e3;
    subplot(411);
    plot(t,abs(rf), 'r-');
    hold on;
    plot(t,real(rf),'b--');
    hold on;
    plot(t,imag(rf), 'g--');
    title('RF Envelope - I/Q');
    ylabel('(Gauss)');
    xlabel('Time [ms]');

    % Gradient
    %
    subplot(412);
    plot(t,g,'b-');
    title('Excitation Gradient');
    ylabel('[G/cm]');
    xlabel('Time [ms]');
    
    switch ptype
        case {'ex', 'se'}
            % abs(Mxy)
            %
            subplot(4,3,7);
            imagesc(f,z,abs(m));
            colormap(gray)
            xlabel('Frequency [Hz]');
            ylabel('Position [cm]');
            title('Magnitude M_{xy}')

            % Angle(Mxy)
            %
            subplot(4,3,8);
            imagesc(f,z,angle(m));
            xlabel('Frequency [Hz]');
            ylabel('Position [cm]');
            title('Phase M_{xy}');
        case {'inv', 'sat'}
            % Mz
            subplot(4,3,[7,8]);
            imagesc(f,z,m);
            colormap(gray)
            xlabel('Frequency [Hz]');
            ylabel('Position [cm]');
            title('M_{z}')
    end    
    
    % Spectral plot at z = 0
    %
    subplot(4,3,9);
    dbw_fine = diff(bw)/499;
    f_fine = [bw(1):dbw_fine:bw(2)];
    m_center = calc_mag(rf_rot,gz_rot+i*gf_rot, 0, f_fine,ptype);
    plot_mag(f_fine, m_center, ptype);
    xlabel('Frequency [Hz]');
    title(sprintf('Spectral Profile - Z = 0'));

    
    % Passband/Stopband plots
    %
    if (nargin < 8) || isempty(f), 
	% Sort frequency from 0 out
	[tmp idx] = sort(abs(f));
	fsort = f(idx);
	
	% Get spectral profile
	m_z0 = calc_mag(rf_rot,gz_rot+i*gf_rot, 0, fsort,ptype);
	
	% Find maximum closest to 0 frequency
	[mval idx] = max(abs(m_z0));
	fplot = fsort(idx);
    end;
    
    nplot = length(fplot);
    for idx = 1:nplot, 
	subplot(4,nplot,3*nplot + idx);
	m_f = calc_mag(rf_rot,gz_rot+i*gf_rot, z, fplot(idx),ptype);
    plot_mag(z,m_f,ptype);
	title(sprintf('%5.1f', fplot(idx)));
	xlabel('Position [cm]');
    end;

end

% Helper functions to avoid lots of "switch" statements

function M = calc_mag(rfs, gs, z, f, ptype)

[a,b] = abr(rfs,gs, z, f);
switch ptype
    case 'ex'
        M = ab2ex(a,b);
    case 'inv'
        M = ab2inv(a,b);
    case 'sat'
        M = ab2sat(a,b);
    case 'se'
        M = ab2se(a, b);
end

end

function plot_mag(x, M, ptype)

switch ptype
    case {'ex', 'se'}
        plot(x,abs(M), 'r-',x,real(M),'b--', x,imag(M), 'g--');
        grid;
        ylabel('M_{xy}');
    case {'inv', 'sat'}
        plot(x,M, 'b-');
        grid;
        ylabel('M_{z}');
end

end