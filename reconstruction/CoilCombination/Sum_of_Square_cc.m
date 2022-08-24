function image_cc = Sum_of_Square_cc(image_mc)
% Sum of Square method for image combination
% Reference:
%   Hall, Emma L., et al. "Methodology for improved detection of low 
%   concentration metabolites in MRS: optimised combination of signals from 
%   multi-element coil arrays." Neuroimage 86 (2014): 35-42.
%
% INPUTS:
%   Image_mc - multi-channel image data. dim:[Nx, Ny, Nz, Ncoil, 1,
%       N_freq, Ndyn]
%       
% OUTPUTS:
%   Image_cc - coil combination Image.
%
% (c) 2017-2018 The Regents of the University of California
% All Rights Reserved.
%
% Written by :Xucheng Zhu
%   
%

image_cc = sqrt(sum(abs(image_mc.^2),4));