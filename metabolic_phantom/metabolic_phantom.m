function [kTRANS, kPL] = metabolic_phantom( nx, ny, nz, kTRANS_low, kTRANS_high, kPL_low, kPL_high )
% METABOLIC_PHANTOM generates standardized 3-dimensional perfusion and metabolism maps for simulated experiments
% 
%   Parameters: 
%       nx          = number of samples in x dimension 
%       ny          = number of samples in y dimension 
%       nz          = number of samples in z dimension 
%       kTRANS_low  = perfusion rate in low perfusion region
%       kTRANS_high = perfusion rate in high perfusion region
%       kPL_low     = metabolic rate in low metabolism regions
%       kPL_high    = metabolic rate in high metabolism regions
% 
%   Outputs: 
%       kTRANS      = generated perfusion map
%       kPL         = generated metabolic map 
% 
%   Author: 
%       John Maidens (maidens@eecs.berkeley.edu)
%       April 2017

    x = linspace(-1, 1, nx);
    y = linspace(-1, 1, ny);
    z = linspace(-1, 1, nz);
    
    [X, Y, Z] = meshgrid(x, y, z);
    
    kTRANS = zeros(size(X)); 
    kTRANS = kTRANS + kTRANS_low*rectangle_shape(X, Y, Z, 0.8, 1.6, 1.6, 0.4, 0, 0);
    kTRANS = kTRANS + kTRANS_high*rectangle_shape(X, Y, Z, 0.8, 1.6, 1.6, -0.4, 0, 0);
    
    kPL = zeros(size(X));
    kPL = kPL + kPL_low*sphere(X, Y, Z, 0.35, 0.45, 0.45, 0) + (kPL_high - kPL_low)*sphere(X, Y, Z, 0.10, 0.45, 0.45, 0);
    kPL = kPL + kPL_low*sphere(X, Y, Z, 0.35, -0.45, 0.45, 0) + (kPL_high - kPL_low)*sphere(X, Y, Z, 0.10, -0.45, 0.45, 0);
    kPL = kPL + kPL_high*sphere(X, Y, Z, 0.35, 0.45, -0.45, 0) + (kPL_low - kPL_high)*sphere(X, Y, Z, 0.10, 0.45, -0.45, 0);
    kPL = kPL + kPL_high*sphere(X, Y, Z, 0.35, -0.45, -0.45, 0) + (kPL_low - kPL_high)*sphere(X, Y, Z, 0.10, -0.45, -0.45, 0);   
    
end


    