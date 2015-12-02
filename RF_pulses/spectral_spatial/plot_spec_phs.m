function plot_spec_phs(f, a, d, w, h, type)
% PLOT_SPEC - Utility to plot frequency specifications
%   
%  plot_spec(f, a, d, type)
%
%  f - frequency band edges
%  a - band amplitudes
%  d - band ripple specs
%  type - line/plotting type (see S in 'help plot')

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

if nargin < 6,
    type = 'k--';
end;
nband = length(f)/2;
subplot(211)
for band = 1:nband,
    idx = [band*2-1 band*2];
    plot(f(idx), abs(a(band))+abs(d(band))*ones(1,2), sprintf('%s',type));
    hold on;
    plot(f(idx), abs(a(band))-abs(d(band))*ones(1,2), sprintf('%s',type));
end;
plot(w,abs(h));

subplot(212)
for band = 1:nband,
    idx = [band*2-1 band*2];
    if abs(a(band)) == 0
        continue;
    end
    plot(f(idx), angle(a(band))+angle(d(band))*ones(1,2), sprintf('%s',type));
    hold on;
    plot(f(idx), angle(a(band))-angle(d(band))*ones(1,2), sprintf('%s',type));
end;
plot(w,angle(h));
