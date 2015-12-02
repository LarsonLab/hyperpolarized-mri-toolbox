 function [d, a, ang] = rf_ripple(de, a, ang, type)
% RF_RIPPLE - Compute polynomial ripple required for various pulse types
%   
% Input
%    - de - vector of ripples required
%    - a - vector of amplitudes in bands
%    - ang - angle of excitation
%    - excitation type: 'ex', 'se', 'sat', 'inv'
%      
% Output
%    - d: ripples for filter design
%    - a: band amplitudes for filter design
%    - ang: angle of excitation

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


nband = length(de);
for band = 1:nband,
    db = de(band);
    ab = a(band);
    switch (type )
        case 'ex'
            % Determine B ripple originating from mxy ripple,
            % find which is most significant, then use this
            % to set the ripple
            %
            mxy_mid = sin(2*asin(sin(ang/2)*ab));
            mxy_pos = max(-1,min(1,db + mxy_mid));
            mxy_neg = max(-1,min(1,-db + mxy_mid));

            B_mid = sin(ang/2) * ab;
            B_pos = sin(asin(mxy_pos)/2);
            B_neg = sin(asin(mxy_neg)/2);

            d(band) = max(abs(B_pos-B_mid), abs(B_neg-B_mid))/sin(ang/2);
        case 'se'
            mxy_mid = (sin(ang/2)*ab)^2;
            mxy_pos = max(0,min(1,db + mxy_mid));
            mxy_neg = max(0,min(1,-db + mxy_mid));

            % check for ripples that will result in B > 1 and scale back
            % a if necessary
            if (mxy_pos >= 1) || (mxy_neg >= 1)
                ab = sqrt(1 - abs(db)) / sin(ang/2);
                mxy_mid = (sin(ang/2)*ab)^2;
                mxy_pos = max(0,min(1,db + mxy_mid));
                mxy_neg = max(0,min(1,-db + mxy_mid));
                a(band) = ab;
            end
            
            B_mid = sin(ang/2) * ab;
            B_pos = sqrt(mxy_pos);
            B_neg = sqrt(mxy_neg);

            d(band) = max(abs(B_pos-B_mid), abs(B_neg-B_mid))/sin(ang/2);
%             if a(band) == 1,
%                 d(band) = de(band)/4;
%             else
%                 d(band) = sqrt(de(band));
%             end;
        case {'inv','sat'}
            mz_mid = 1 - 2* (sin(ang/2)*ab)^2;
            mz_pos = max(-1,min(1,db + mz_mid));
            mz_neg = max(-1,min(1,-db + mz_mid));

            % check for ripples that will result in B > 1 and scale back
            % a if necessary
            if (mz_neg <= -1) || (mz_pos <= -1)
                ab = sqrt(1 - abs(db)/2) / sin(ang/2);
                mz_mid = 1 - 2* (sin(ang/2)*ab)^2;
                mz_pos = max(-1,min(1,db + mz_mid));
                mz_neg = max(-1,min(1,-db + mz_mid));
                a(band) = ab;
            end
            
            B_mid = sin(ang/2) * ab;
            B_pos = sin(acos(mz_pos)/2);
            B_neg = sin(acos(mz_neg)/2);
%            B_pos = sqrt((1-mz_pos)/2);  % same as above
%            B_neg = sqrt((1-mz_neg)/2);

            d(band) = max(abs(B_pos-B_mid), abs(B_neg-B_mid))/sin(ang/2);
%             if a(band) == 1,
%                 d(band) = de(band)/8;
%             else
%                 d(band) = sqrt(de(band)/2);
%             end;
%         case 'sat'
%             if a(band) == 1,
%                 d(band) = de(band)/2;
%             else
%                 d(band) = sqrt(de(band));
%             end;

    end;
end;
    

% scale a if the peak angle has been reduced
if max(a) ~= 1
    max_a = max(a);
    a = a / max_a;
    ang = 2*asin(sin(ang/2)*max_a);
    d = d / max_a;
end