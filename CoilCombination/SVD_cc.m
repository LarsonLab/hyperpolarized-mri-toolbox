function Image_cc = SVD_cc(Image_mc)
% SVD based Coil Combination for Hyperpolarized MRI data
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

[Nx, Ny, Nz, Ncoil, Nmap, N_freq, Ndyn] = size(Image_mc);
Image_mc = prewhitening_cc(Image_mc);
Image_cc = zeros([Nx, Ny, Nz, 1, Nmap, N_freq, Ndyn]);

for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            for t = 1:Ndyn
                tmp_spectrum = Image_mc(i,j,k,:,1,:,t);
                tmp_spectrum = squeeze(tmp_spectrum);
                [U,S,V] = svd(tmp_spectrum,'econ');
                S = wthresh(S,1);
                Image_cc(i,j,k,1,1,:,t) = sum(U*S*V',2);
            end
        end
    end
end
