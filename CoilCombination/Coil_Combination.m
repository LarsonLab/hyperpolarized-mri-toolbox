function [Image_cc, Smap] = Coil_Combination(Image_mc, cc_method_flag, prewhitening_flag)
% Coil Combination for Hyperpolarized MRI data
% Reference:
%   Hall, Emma L., et al. "Methodology for improved detection of low 
%   concentration metabolites in MRS: optimised combination of signals from 
%   multi-element coil arrays." Neuroimage 86 (2014): 35-42.
%
% INPUTS:
%   Image_mc - multi-channel image data. dim:[Nx, Ny, Nz, Ncoil, 1,
%       N_freq, Ndyn]
%   cc_method_flag -
%       0 - sum of square (SOS).
%       1 - Singular value decomposition(SVD). Only works for spetroscopy
%       data
%       2 - Reference Peak method.
%   prewhitening_flag - 
%       0 - No prewhitening
%       1 - Prewhitening based on edge of image
%       
% OUTPUTS:
%   Image_cc - coil combination Image.
%   Smap - sensitivity maps, empty for 0,1,2 results.
%
% (c) 2017-2018 The Regents of the University of California
% All Rights Reserved.
%
% Written by :Xucheng Zhu
%   
%

if nargin < 2
    cc_method_flag = 0;
    prewhitening_flag = 0;
end

if nargin < 3
    prewhitening_flag = 0;
end

if prewhitening_flag
    % assuming the last time point is noise
    [Image_cc, Correlation_Matrix] = prewhitening_cc(Image_cc);
end

[Nx, Ny, Nz, Ncoil, Nmap, N_freq, Ndyn] = size(Image_mc);

switch cc_method_flag
    case 0
        Image_cc = Sum_of_Square_cc(Image_mc);
        Smap = [];
    case 1
        Image_cc = SVD_cc(Image_mc);
        Smap = [];
    case 2
        [Image_cc, Smap] = RefPeak_cc(Image_mc);
end