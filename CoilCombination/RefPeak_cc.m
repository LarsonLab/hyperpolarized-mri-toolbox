function [Image_cc,Smap] = RefPeak_cc(Image_mc)
% Reference Peak method for image combination
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
%   Smap - sensitivity map
% (c) 2017-2018 The Regents of the University of California
% All Rights Reserved.
%
% Written by :Xucheng Zhu
%   
%

[Nx, Ny, Nz, Ncoil, Nmap, N_freq, Ndyn] = size(Image_mc);
% sum up data
Image_t = sum(Image_mc,7);
Spectrum_t = squeeze(sum(sum(sum(sum(Image_t,1),2),3),4));
[~,peak_pos] = max(squeeze(abs(Spectrum_t)));
image_ref = Image_t(:,:,:,:,1,peak_pos);

Image_sos = sqrt(abs(sum(image_ref.*conj(image_ref),4)));
Smap = Image_t./repmat(Image_sos,1,1,1,Ncoil,1,1,1);

Image_cc = sum(Image_mc.*repmat(conj(Smap),1,1,1,1,1,N_freq,Ndyn),4);