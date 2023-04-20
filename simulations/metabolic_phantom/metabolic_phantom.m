function [kTRANS, kPL] = metabolic_phantom( nx, ny, nz, kTRANS_low, kTRANS_high, kPL_low, kPL_high, linear_kTRANS_gradient, linear_kPL_gradient)
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
%       linear_kTRANS_gradient = boolean flag for whether kTRANS should be set to a gradient
%       linear_kPL_gradient    = boolean flag for whether kPL should be set to a gradient
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
    
    % option for 2D phantom
    if nz > 1
        z = linspace(-1, 1, nz);
    else
        z = 0;
    end
    [X, Y, Z] = meshgrid(x, y, z);
    
    kTRANS = zeros(size(X)); 
    
    if (linear_kTRANS_gradient)
        scale = 0.5*(kTRANS_low - kTRANS_high)*X + 0.5*(kTRANS_low + kTRANS_high);
        kTRANS = kTRANS + scale.*rectangle_shape(X, Y, Z, 1.6, 1.6, 1.6, 0, 0, 0);
    else
        kTRANS = kTRANS + kTRANS_low*rectangle_shape(X, Y, Z, 0.8, 1.6, 1.6, 0.4, 0, 0);
        kTRANS = kTRANS + kTRANS_high*rectangle_shape(X, Y, Z, 0.8, 1.6, 1.6, -0.4, 0, 0);
    end
    
    kPL = zeros(size(X));
    if (linear_kPL_gradient)
        scale = 0.5*(kPL_low - kPL_high)*Y + 0.5*(kPL_low + kPL_high);
        kPL = kPL + scale.*sphere_fcn(X, Y, Z, 0.35, 0.45, 0.45, 0) + (kPL_high - scale).*sphere_fcn(X, Y, Z, 0.10, 0.45, 0.45, 0);
        kPL = kPL + scale.*sphere_fcn(X, Y, Z, 0.35, -0.45, 0.45, 0) + (kPL_high - scale).*sphere_fcn(X, Y, Z, 0.10, -0.45, 0.45, 0);
        kPL = kPL + scale.*sphere_fcn(X, Y, Z, 0.35, 0.45, -0.45, 0) + (kPL_low - scale).*sphere_fcn(X, Y, Z, 0.10, 0.45, -0.45, 0);
        kPL = kPL + scale.*sphere_fcn(X, Y, Z, 0.35, -0.45, -0.45, 0) + (kPL_low - scale).*sphere_fcn(X, Y, Z, 0.10, -0.45, -0.45, 0);   
    else
        kPL = kPL + kPL_low*sphere_fcn(X, Y, Z, 0.35, 0.45, 0.45, 0) + (kPL_high - kPL_low)*sphere_fcn(X, Y, Z, 0.10, 0.45, 0.45, 0);
        kPL = kPL + kPL_low*sphere_fcn(X, Y, Z, 0.35, -0.45, 0.45, 0) + (kPL_high - kPL_low)*sphere_fcn(X, Y, Z, 0.10, -0.45, 0.45, 0);
        kPL = kPL + kPL_high*sphere_fcn(X, Y, Z, 0.35, 0.45, -0.45, 0) + (kPL_low - kPL_high)*sphere_fcn(X, Y, Z, 0.10, 0.45, -0.45, 0);
        kPL = kPL + kPL_high*sphere_fcn(X, Y, Z, 0.35, -0.45, -0.45, 0) + (kPL_low - kPL_high)*sphere_fcn(X, Y, Z, 0.10, -0.45, -0.45, 0);   
    end
end


    