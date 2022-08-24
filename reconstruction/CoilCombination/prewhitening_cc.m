function [Image_w, Correlation_Matrix] = prewhitening_cc(Image_mc, Noise_data)
% Prewhitening step for Image_data
% Reference:
%   Hall, Emma L., et al. "Methodology for improved detection of low 
%   concentration metabolites in MRS: optimised combination of signals from 
%   multi-element coil arrays." Neuroimage 86 (2014): 35-42.
%
% INPUTS:
%   Image_mc - multi-channel image data. dim:[Nx, Ny, Nz, Ncoil, 1,
%       N_freq, Ndyn]
%
%   Noise_data - Optional argument to specify noise-only data. If not
%   supplied, will default to using the final timepoint from all
%   frequencies and all slices from Image_mc.
%       
% OUTPUTS:
%   Image_w - Prewhitened image data.
%   Correlation_Matrix - Correlation
%
% (c) 2017-2018 The Regents of the University of California
% All Rights Reserved.
%
% Written by :Xucheng Zhu
%   
%
[Nx, Ny, Nz, Ncoil, Nmap, N_freq, Ndyn] = size(Image_mc);
if nargin == 1
    Noise = Image_mc(1,:,:,:,1,:,end);
else
    Noise = Noise_data;
end
Noise = permute(Noise,[4,1,2,3,5,6,7]);
Nsample = numel(Noise)/Ncoil;
% TODO: add minimum points for prewhitening calibration
Noise = reshape(Noise,Ncoil,Nsample);

Correlation_Matrix = abs(Noise*Noise')/(Nsample+eps);

Sigma = permute(diag(sqrt(Correlation_Matrix)),[2 3 4 1]);
Sigma = repmat(Sigma,Nx, Ny, Nz, 1, Nmap, N_freq, Ndyn);

Image_w = Image_mc./(Sigma+eps);



